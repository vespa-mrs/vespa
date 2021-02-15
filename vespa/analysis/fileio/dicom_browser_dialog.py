# Python modules
from __future__ import division, print_function
import os
import os.path

# 3rd party modules
import wx
import pydicom
import numpy as np

from wx.lib.embeddedimage import PyEmbeddedImage


__version__ = '1.4.0'


# These content type constants are arbitrary and may change.
# However bool(NONE) is guaranteed to be False
TYPE_NONE = 0
TYPE_IMAGE = 1
TYPE_SPECTROSCOPY = 2




# -----------------------------------------------------------------------------
# The GUI code - auto-generated in wxGlade and copied into this section

class TheDialog(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: TheDialog.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, *args, **kwds)
        self.SetSize((653, 758))
        self.TextImportDirectory = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.buttonChooseDirectory = wx.Button(self, wx.ID_ANY, "Browse ...")
        self.SplitterFilesTags = wx.SplitterWindow(self, wx.ID_ANY)
        self.PanelFiles = wx.ScrolledWindow(self.SplitterFilesTags, wx.ID_ANY, style=wx.TAB_TRAVERSAL)
        self.TreeFiles = wx.TreeCtrl(self.PanelFiles, wx.ID_ANY, style=wx.BORDER_SIMPLE | wx.TR_DEFAULT_STYLE | wx.TR_HAS_BUTTONS | wx.TR_HIDE_ROOT | wx.TR_LINES_AT_ROOT | wx.TR_NO_LINES)
        self.PanelTags = wx.ScrolledWindow(self.SplitterFilesTags, wx.ID_ANY, style=wx.TAB_TRAVERSAL)
        self.TreeTags = wx.TreeCtrl(self.PanelTags, wx.ID_ANY)
        self.labelDeveloperComment = wx.StaticText(self, wx.ID_ANY, "The preview bitmap and the Cancel and Open \nbuttons are added in dialog __init__(), not here.")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.on_choose_directory, self.buttonChooseDirectory)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.on_tree_item_activated, self.TreeFiles)
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.on_tree_selection_changed, self.TreeFiles)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: TheDialog.__set_properties
        self.SetTitle("dialog_1")
        self.SetSize((653, 758))
        self.buttonChooseDirectory.SetFocus()
        self.PanelFiles.SetScrollRate(10, 10)
        self.PanelTags.SetScrollRate(10, 10)
        self.SplitterFilesTags.SetMinimumPaneSize(20)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: TheDialog.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.FlexGridSizer(1, 4, 0, 0)
        sizer_3 = wx.StaticBoxSizer(wx.StaticBox(self, wx.ID_ANY, "DICOM Import Directory"), wx.VERTICAL)
        sizer_4 = wx.StaticBoxSizer(wx.StaticBox(self.PanelTags, wx.ID_ANY, "DICOM Tags"), wx.VERTICAL)
        sizer_2 = wx.StaticBoxSizer(wx.StaticBox(self.PanelFiles, wx.ID_ANY, "Subjects"), wx.VERTICAL)
        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5.Add(self.TextImportDirectory, 1, wx.ALIGN_CENTER_VERTICAL | wx.TOP, 1)
        sizer_5.Add(self.buttonChooseDirectory, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_3.Add(sizer_5, 0, wx.ALL | wx.EXPAND, 0)
        sizer_2.Add(self.TreeFiles, 5, wx.ALL | wx.EXPAND, 2)
        self.PanelFiles.SetSizer(sizer_2)
        sizer_4.Add(self.TreeTags, 1, wx.EXPAND, 0)
        self.PanelTags.SetSizer(sizer_4)
        self.SplitterFilesTags.SplitVertically(self.PanelFiles, self.PanelTags)
        sizer_3.Add(self.SplitterFilesTags, 1, wx.EXPAND, 0)
        sizer_1.Add(sizer_3, 1, wx.EXPAND | wx.LEFT | wx.TOP, 2)
        grid_sizer_1.Add(self.labelDeveloperComment, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.AddGrowableRow(0)
        grid_sizer_1.AddGrowableCol(0)
        grid_sizer_1.AddGrowableCol(1)
        sizer_1.Add(grid_sizer_1, 0, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

    def on_choose_directory(self, event):  # wxGlade: TheDialog.<event_handler>
        print("Event handler 'on_choose_directory' not implemented!")
        event.Skip()

    def on_tree_item_activated(self, event):  # wxGlade: TheDialog.<event_handler>
        print("Event handler 'on_tree_item_activated' not implemented!")
        event.Skip()

    def on_tree_selection_changed(self, event):  # wxGlade: TheDialog.<event_handler>
        print("Event handler 'on_tree_selection_changed' not implemented!")
        event.Skip()


# -----------------------------------------------------------------------------
# The actual dialog widget

class DialogDicomBrowser(TheDialog):
    """
    Opens a dialog that allows the user to select a directory containing DICOM
    files. Displays the files in a tree and allows the user to select one or
    multiple files.

    A Series folder can be selected to select all files in that series, and a
    study folder can be selected to return all series/files under it.

    Multiple series can be selected at a time and the User can mix and match
    individual files and series folders so long as the 'multi_select' option
    is True.

    After the dialog closes, the selected filenames are stored in the 'items'
    attribute of this class as either 1) the string name of each selected file,
    or 2) the DicomFileInfo object(s) for each selected file. Each
    DicomFileInfo contains a reference to the actual pydicom.Dataset object
    for each file. Note that for the sake of consistency, dialog.items is a
    list even if the dialog is invoked with the multi_select option set to
    False.

    If the user doesn't select a file, dialog.items is an empty list.

    Options:

    A preview image can be displayed in lower left if image element selected
    is an image file that our code understands. If more than one image element
    is selected, no preview image is shown. Use the 'preview_size' keyword to
    control this functionality

    Middle panel can be split into a second (right hand side) TreeCtrl that
    displays the DICOM tags for the image element selected in the left TreeCtrl.
    If more than one image element is selected, no tags are displayed. Use the
    'show_tags' keyword to control this functionality.

    """
    content_types = { TYPE_NONE          : "",
                      TYPE_IMAGE         : "Image",
                      TYPE_SPECTROSCOPY  : "Spectroscopy",
                    }

    def __init__(self, parent, preview_size=(64, 64),
                               show_tags=None,
                               multi_select=False,
                               return_objects=False,
                               default_path=''):
        """
        Creates an instance of the dialog.

        Args:
            parent (wxPythonID): parent GUI element to which the dialog attaches.

            preview_size (None, 2-tuple): default=(64,64). If preview_size is
                None, the dialog won't display previews of image files.
                Otherwise the preview size should be a tuple that describes the
                preview size.

            show_tags (None, bool): default None. If set to True/False a 'Show
                tags' checkbox is displayed in the bottom row and set to the
                given value. When set True, the middle panel is split into a
                second TreeCtrl on the right hand side to show the DICOM tags
                for the image element selected on the left side tree. Set to
                None to not show the 'Show tags' checkbox and hide the DICOM
                tags tree.

            multi_select (bool): default False. If False the tree only allows
                single selection. Note that this flag must be True to allow
                users to select a Series folder since this will return multiple
                files (depending on what is in the series).

            return_objects (bool): default False. If False, only the file names
                of the dicom files selected are returned. If True, then the
                pydicom dataset references for each file are returned.

            default_path (string): default ''. If set, dialog will initialize
                with any results from the path, if it exists and there are
                DICOM files there.

        Returns:
            After the dialog closes, the selected filenames are stored in the
            'items' attribute of this class as either 1) the string name of
            each selected file, or 2) the DicomFileInfo object(s) for each
            selected file. Each DicomFileInfo contains a reference to the
            actual pydicom.Dataset object for each file. Note that for the
            sake of consistency, dialog.items is a list even if the dialog
            is invoked with the multi_select option set to False.

        """
        style = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        TheDialog.__init__(self, parent, -1, style=style)

        # Note. windows appears to need 16x16 size or does not appear
        self.SetTitle("DICOM Browser")
        self.SetIcon(top_icon.GetIcon())

        self.multi_select = multi_select
        self.return_objects = return_objects

        if multi_select:
            style = self.TreeFiles.GetWindowStyle()
            style &= ~wx.TR_SINGLE
            style |= wx.TR_MULTIPLE
            self.TreeFiles.SetWindowStyle(style)

        # I added a comment to the dialog for wxGlade users to explain the
        # lack of buttons. It also serves a purpose here, which is that is
        # allows me to grab its containing (grid) sizer which I need here.
        # I don't want the comment shown at runtime, however, so I delete it.
        grid_sizer = self.labelDeveloperComment.GetContainingSizer()
        self.labelDeveloperComment.Destroy()

        # Add the preview bitmap as first item in the grid sizer.
        if preview_size:
            self.preview_size = preview_size
            self.preview = wx.StaticBitmap(self, -1, wx.NullBitmap, size=preview_size)
            grid_sizer.Add(self.preview, 0, wx.LEFT | wx.BOTTOM, 10)
        else:
            grid_sizer.Add((20,20), 0, 0, 0)    # no preview, but need for proper spacing
            self.preview = None

        self.LabelCount = wx.StaticText(self, wx.ID_ANY, "file count goes here...")
        grid_sizer.Add(self.LabelCount, 0, wx.BOTTOM | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 14)

        self.show_tags = show_tags
        if show_tags is not None:
            self.CheckShowTags = wx.CheckBox(self, wx.ID_ANY, " Show Tags ")
            grid_sizer.Add(self.CheckShowTags, 0, wx.BOTTOM | wx.ALIGN_BOTTOM | wx.ALIGN_RIGHT, 14)
            self.Bind(wx.EVT_CHECKBOX, self.on_show_tags, self.CheckShowTags)
            if show_tags:
                self.SplitterFilesTags.Unsplit()

        # The Open and Cancel buttons appear in a different order on Windows
        # than on OS X & GTK, but the wx.StdDialogButtonSizer handles that
        # for me.
        # Also, using wx.ID_OPEN and wx.ID_CANCEL ensures that these buttons
        # get the appropriate artwork under GTK.
        # Last but not least, a button with id == ID_CANCEL is automatically
        # bound to self.Close()
        button_sizer = wx.StdDialogButtonSizer()
        self.buttonOpen = wx.Button(self, wx.ID_OPEN, "Open")
        self.buttonCancel = wx.Button(self, wx.ID_CANCEL, "Cancel")
        button_sizer.AddButton(self.buttonOpen)
        button_sizer.AddButton(self.buttonCancel)
        button_sizer.SetAffirmativeButton(self.buttonOpen)
        self.Bind(wx.EVT_BUTTON, self.on_open, self.buttonOpen)
        self.buttonOpen.SetDefault()
        # Open button enabled only when a valid item is selected.
        self.buttonOpen.Enable(False)
        
        button_sizer.Realize()
        grid_sizer.Add(button_sizer, 0, wx.BOTTOM|wx.ALIGN_BOTTOM, 10)

        self.Layout()

        # I add a fake root to the tree. The fake root will be invisible.
        # Allows me to display multiple sets of DICOM data as siblings.
        self.TreeFiles.AddRoot("")

        self._user_home = wx.GetUserHome()

        self.items = []
        self.path = default_path
        self.tags_root = False
        self.sash_pos = 300
        self._eat_a_selection_event = False

        self._set_status('')

        # Here I set up the icons used in the tree.
        image_size = (16, 16)
        images = wx.ImageList(image_size[0], image_size[1])

        icon1 = wx.ArtProvider.GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, image_size)
        icon2 = wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, image_size)
        self.icon_index_user   = images.Add(user_icon.GetIcon())
        self.icon_index_book   = images.Add(book_icon.GetIcon())
        self.icon_index_folder = images.Add(icon1)
        self.icon_index_file   = images.Add(icon2)

        self.TreeFiles.AssignImageList(images)   # order matters?

        self.Center()

        # here we set up the Tags display and populate the GUI as appropriate
        if self.show_tags is not None:
            self.CheckShowTags.SetValue(self.show_tags)
            if self.show_tags:
                self.SplitterFilesTags.SplitVertically(self.PanelFiles, self.PanelTags)
                self.sash_pos = self.SplitterFilesTags.GetSashPosition()
            else:
                self.SplitterFilesTags.Unsplit(self.PanelTags)
        else:
            self.SplitterFilesTags.Unsplit(self.PanelTags)

        if self.path:
            self.TextImportDirectory.SetValue(self.path)
            self._populate_files_tree()

    #--------------------------------------------------------------------------
    # control settings for the Dialog's Frame

    __doc = """User's currently selected path which is reflected in the dialog's title bar."""
    def __get_path(self):
        return self._path
    def __set_path(self, path):
        self._path = path
        if not path:
            path = "[no directory selected]"
        elif path.startswith(self._user_home):
            # Trim off the home dir + separator to save space
            path = path[len(self._user_home) + 1:]
            
        # bjs - commented this out in this version since I'm using a widget
        #       TextCtrl to display the selected directory
        #self.SetTitle("DICOM Browser - " + path)
    path = property(__get_path, __set_path, doc=__doc)


    @property
    def filenames(self):
        if self.items:
            if hasattr(self.items[0], "filename"):
                return [item.filename for item in self.items]
            else:
                return self.items
        else:
            return []


    #--------------------------------------------------------------------------
    # Overloaded methods

    def file_filter(self, dicom_file):
        """
        Called for each file to check if should be included in the tree.
        If the function returns False, the file is excluded.

        This version of the function excludes no files because always True.
        Subclass this dialog and override this function for custom filter.
        """
        return True


    #--------------------------------------------------------------------------
    # Event Handlers

    def on_tree_selection_changing(self, event):
        if self._eat_a_selection_event:
            # Eat this event, but not subsequent ones.
            # See on_tree_selection_changed().
            event.Veto()
            self._eat_a_selection_event = False
        else:
            event.Skip()

    def on_tree_selection_changed(self, event):

        item_ids = self._get_selections()

        if self.multi_select and (len(item_ids) == 1):
            # Multi-select is enabled, and a single item was selected. If it's
            # a header (a non-leaf node), select all children.
            item_id = item_ids[0]
            if self.TreeFiles.ItemHasChildren(item_id):
                # Select all kids. Selecting items tends to scroll the list. We
                # save the first item so we can restore the current scroll
                # position when selecting is done. We also disable selection
                # notifications while doing this -- it speeds things up
                # considerably.
                first_visible = self.TreeFiles.GetFirstVisibleItem()
                self.Unbind(wx.EVT_TREE_SEL_CHANGED, self.TreeFiles)
                self._select_offspring(item_id)
                self.Bind(wx.EVT_TREE_SEL_CHANGED, self.on_tree_selection_changed, self.TreeFiles)

                # Strange but true: under OS X & GTK, an EVT_LEFT_UP event
                # occurs after this event, and that triggers an item selection
                # which undoes the work of _select_offspring(). To work around
                # this, I set a flag which eats the next selection event.
                if wx.Platform != "__WXMSW__":
                    self._eat_a_selection_event = True

                self.TreeFiles.ScrollTo(first_visible)
            # else:
            #     If only single-selection is enabled, or multiple items have been
            #     selected, don't enable the header-implies-selection behavior.

        filenames = self._get_selected_filenames()
        if not filenames: return

        self.buttonOpen.Enable(bool(filenames))  # enabled while 1+ filename selected

        self._set_status("{0} of {1} files selected.".format(len(filenames), self._nfiles))


        if self.show_tags is not None:
            if self.show_tags:
                self.update_tags()

        if self.preview:
            bitmap = wx.NullBitmap

            # Preview gets pretty confusing in multi-select mode. The code
            # below handles both single- and multi-select.

            # Preview is ambiguous when multiple items are selected, so we
            # turn it off in that case. Note that with creative clicking
            # it is possible to have > 1 items selected even when
            # len(self._get_selected_filenames()) == 1. Therefore it's not
            # sufficient for me to check len(filenames) here, I have to check
            # how many items are selected.

            # In addition, note that this event fires when an item is
            # selected or *de*selected; in the latter case the item returned
            # by event.GetItem() won't be selected. Imagine if the user
            # selects items 5 - 8, and then deselects 7, 6, and 5 in that
            # order. When #5 is deselected, there will only be one item
            # selected (#8) and so its preview should appear. However,
            # event.GetItem() will point to item #5. Be careful!

            item_ids = self._get_selections()
            if len(item_ids) == 1:
                # OK, only one item is selected. Now ensure that it has
                # a filename associated.
                items = self.TreeFiles.GetItemData(item_ids[0])
                if items:
                    if len(items) == 1:
                        for item in items:
                            filename = item.filename
                            try:
                                image = item.get_image()
                            except NotImplementedError:
                                image = None                # image in a unknown format

                            if image:
                                size = self.preview.GetSize()
                                image.Rescale(size[0], size[1])
                                bitmap = wx.Bitmap(image)
                # else:
                # The item selected is not a leaf node
            # else:
            # There's more than one item selected

            self.preview.SetBitmap(bitmap)


    def on_choose_directory(self, event):
        self.path = wx.DirSelector("Select a directory of DICOM files",
                                   self._suggest_a_start_path())
        if self.path:
            self.TextImportDirectory.SetValue(self.path)
            if self.tags_root:
                self.TreeTags.DeleteChildren(self.tags_root)
            self._populate_files_tree()


    def on_show_tags(self, event):
        self.show_tags = event.GetEventObject().GetValue()
        if self.show_tags:
            self.SplitterFilesTags.SplitVertically(self.PanelFiles,self.PanelTags)
            self.SplitterFilesTags.SetSashPosition(self.sash_pos, redraw=True)
            self.update_tags()
        else:
            self.sash_pos = self.SplitterFilesTags.GetSashPosition()
            self.SplitterFilesTags.Unsplit(self.PanelTags)


    def on_tree_item_activated(self, event):
        """
        'Activated' == a double click,  I don't want a double click to do
         anything if the currently selected item isn't a leaf node.

        """
        item_id = event.GetItem()
        if item_id and (not self.TreeFiles.ItemHasChildren(item_id)):
            self.on_open(event)


    def on_open(self, event):
        item_ids = self._get_selections()

        self.items = []
        if item_ids:
            for item_id in item_ids:
                if not self.TreeFiles.ItemHasChildren(item_id):     # exclude Series/Study leaves
                    items = self.TreeFiles.GetItemData(item_id)
                    if items:
                        for item in items:
                            if self.return_objects:
                                self.items.append(item)
                            else:
                                self.items.append(item.filename)

        if not self.multi_select:       # final check for multi_select consistency
            if len(self.items) > 1:
                self.items = []

        self.Close()


    def update_tags(self):
        """generic call that can be used in multiple locations"""
        item_ids = self._get_selections()
        if len(item_ids) == 1:
            # OK, only 1 item is selected. Ensure it has a dataset associated.
            items = self.TreeFiles.GetItemData(item_ids[-1])
            if items:
                if len(items) == 1:
                    for item in items:
                        ds = item.dataset
                        self._populate_tags_tree(ds)
        else:
            if self.tags_root:
                self.TreeTags.DeleteChildren(self.tags_root)


    #--------------------------------------------------------------------------
    # Public and Utility methods

    def get_item_description(self, dicom_file):
        """
        Builds the description string for a leaf item in the list and
        returns that string along with an index pointing to the appropriate
        icon for this item in the tree's ImageList.

        If you want a different description, subclass this dialog and
        override this function. Return -1 for the icon index if you don't
        want an icon.
        """
        description = "%d %s" % (dicom_file.instance_number,
                                 os.path.split(dicom_file.filename)[1])

        # Siemens-specific code -- if there's a content type, display
        # that in the item description
        if hasattr(dicom_file, "content_type"):
            content_type = DialogDicomBrowser.content_types[dicom_file.content_type]
            description = "[%s] %s" % (content_type, description)

        return description, self.icon_index_file


    def is_dicom(self, filename):
        """
        Returns True if the file in question is a DICOM file, else False.

        - a DICOM file starts with 128 reserved bytes followed by "DICM".
        - ref: DICOM spec, Part 10: Media Storage and File Format for Media
               Interchange, 7.1 DICOM FILE META INFORMATION

        """
        if os.path.isfile(filename):
            f = open(filename, "rb")
            s = f.read(132)
            if isinstance(s, (bytes, bytearray)):
                try:
                    s = s.decode('utf-8')
                except:
                    # bugfix - was trying to read PNG and found a char utf8 did not like
                    try:
                        s = s.decode('utf-16')
                    except:
                        return False
            f.close()
            return s.endswith("DICM")
        else:
            return False


    def get_all_files(self, path):
        """
        Returns list of DicomFileInfo instances for all DICOM files in path.
        To sort the files, call the sort() method of the returned list.

        """
        files = []
        for f in self.get_files(path):
            files.append(f)
        return files


    def get_files(self, path):
        """
        Returns the DICOM files in the directory indicated by the path.

        This is a Python generator for DicomFileInfo instances (or a subclass
        thereof) that are returned one by one. Use it in a loop like so:

            for dicom_file in get_files(the_path):
                do_something(dicom_file)

        See also get_all_files().

        """
        # get directory contents, fully-qualified filenames, remove non-files
        filenames = os.listdir(path)
        filenames = [os.path.join(path, filename) for filename in filenames]
        filenames = [filename for filename in filenames if os.path.isfile(filename)]

        for filename in filenames:
            if self.is_dicom(filename):
                dataset = pydicom.read_file(filename)
                dataset.decode()  # change strings to unicode

                if ("Manufacturer" in dataset) and \
                   ("SIEMENS" in dataset.Manufacturer.upper()):
                    df = DicomFileSiemens(dataset, filename)
                else:
                    df = DicomFileInfo(dataset, filename)
                yield df
            # else:
                # Not a DICOM file; ignore it.


    #--------------------------------------------------------------------------
    # Internal methods

    def _suggest_a_start_path(self):
        """Returns a path for the user to begin browsing for DICOM files."""
        # If the user already selected a path, I start there.
        path = self.path

        # Note that wx.StandardPaths may return paths that don't exist!
        if not os.path.exists(path):
            path = wx.StandardPaths.Get().GetDocumentsDir()

        if not os.path.exists(path):
            path = wx.StandardPaths.Get().GetHomeDir()

        if not os.path.exists(path):
            path = os.getcwd()

        return path


    def _has_series_selected(self):
        """ Checks to see if any of the tree items selected have children """
        has_series = False
        item_ids = self.TreeFiles.GetSelections()
        if item_ids:
            for item_id in item_ids:
                if item_id.IsOk():
                    if self.TreeFiles.ItemHasChildren(item_id):
                        has_series = True
        return has_series


    def _get_selections(self):
        """ Always returns the tree item(s) selected in a list """
        if self.multi_select:
            item_ids = self.TreeFiles.GetSelections()
        else:
            item_ids = [self.TreeFiles.GetSelection(),]

        return item_ids


    def _get_selected_filenames(self):
        """
        This method gets the filenames stored in the user data slot of
        the tree items selected. 
        
        It also obeys the "multi_select" flag. So if multi_select == False 
        and more than one filename is reported from the selections, we null 
        the filenames array to indicate that an inappropriate node selection 
        has been made given the flag status.
        
        """
        item_ids = self._get_selections()

        filenames = []
        for item_id in item_ids:
            if not self.TreeFiles.ItemHasChildren(item_id):
                items = self.TreeFiles.GetItemData(item_id)
                if items:
                    for item in items:
                        filenames.append(item.filename)

        # return no filenames from Series node if multi_select is False
        if not self.multi_select and len(filenames) > 1:
            filenames = []

        return filenames


    def _select_offspring(self, item_id):
        """Given a tree item id, recursively selects that item's children,
        including grandchildren, great-grandchildren, etc.
        """
        kid_id, cookie = self.TreeFiles.GetFirstChild(item_id)
        while kid_id.IsOk():
            self.TreeFiles.SelectItem(kid_id)
            if self.TreeFiles.ItemHasChildren(kid_id):
                self._select_offspring(kid_id)
            kid_id, cookie = self.TreeFiles.GetNextChild(item_id, cookie)


    def _populate_files_tree(self):
        """ Populate the TreeFiles element from selected directory """
        wx.SetCursor(wx.HOURGLASS_CURSOR)

        # Reset controls
        self.TreeFiles.DeleteChildren(self.TreeFiles.GetRootItem())
        self.buttonOpen.Enable(False)
        self._set_status("Scanning for DICOM files...")
        self.items = []
        if self.preview:
            self.preview.SetBitmap(wx.NullBitmap)

        # Under Windows & GTK, this dialog needs a repaint after the file
        # open dialog disappears. Unless I force one manually here, it won't
        # happen until this lengthy function completes. It doesn't always
        # get the job done under GTK. :-/
        self.Update()

        # Per issue 58, pydicom doesn't handle paths that are Unicode strings
        # properly (even if the contents are 100% ASCII), so I make sure the
        # path I pass is a string.
        # http://code.google.com/p/pydicom/issues/detail?id=58&can=1&q=defer_size&colspec=ID%20Type%20Status%20Priority%20Milestone%20Owner%20Summary%20Difficulty
        files = []
        self._nfiles = 0
        for f in self.get_files(str(self.path)):
            if self.file_filter(f):
                files.append(f)
                self._nfiles += 1
                self._set_status("Found %d DICOM files." % self._nfiles)

        self._set_status("Found %d DICOM files." % self._nfiles)

        files.sort()

        # Here I iterate through my list and represent it hierachically.
        patient_id = None
        patient_name = None
        study_id = None
        study_description = None
        series_number = None
        series_description = None
        root = self.TreeFiles.GetRootItem()
        current_series = None
        current_series_files = None
        for df in files:
            if df.patient_id != patient_id:
                # New patient
                study_id = None
                series_number = None
                patient_id = df.patient_id
                patient_name = df.patient_name
                s = "%s (%s)" % (df.patient_name, df.patient_id)
                current_patient = self.TreeFiles.AppendItem(root, s, self.icon_index_user)

            if df.study_id != study_id:
                # New study
                series_number = None
                study_id = df.study_id
                study_description = df.study_description
                s = "Study %s : %s" % (df.study_id, df.study_description )
                current_study = self.TreeFiles.AppendItem(current_patient, s, self.icon_index_book)

            if df.series_number != series_number:
                # set list of all files for last series into the last series 
                # tree element before creating a new series node
                if current_series_files:
                    self.TreeFiles.SetItemData(current_series, current_series_files)
                # start a new series node
                series_number = df.series_number
                series_description = df.series_description
                s = "Series %d : %s" % (df.series_number, df.series_description)
                current_series = self.TreeFiles.AppendItem(current_study, s, self.icon_index_folder)
                current_series_files = []

            description, icon_index = self.get_item_description(df)

            item_id = self.TreeFiles.AppendItem(current_series, description, icon_index)
            # Each leaf node's data is set to the associated DicomFileInfo instance.
            self.TreeFiles.SetItemData(item_id, [df])
            current_series_files.append(df)

        # After the last series is appended to tree there is no new series 
        # number to trigger loading the list to the series node, so we 
        # force it to happen now.
        if current_series is not None:
            self.TreeFiles.SetItemData(current_series, current_series_files)

        # could do this, but annoying for large directories
        # self.TreeFiles.ExpandAll()

        if files:
            # Under Windows the list is scrolled to the bottom, I think
            # because of the call to ExpandAll() above. We force the first
            # item to be visible.
            first_child, _ = self.TreeFiles.GetFirstChild(root)
            self.TreeFiles.ScrollTo(first_child)
        else:
            self._set_status("No DICOM files found.")

        wx.SetCursor(wx.NullCursor)

    def _populate_tags_tree(self, ds):
        """ Populate TreeTags element with the [desired] dataset values"""
        self.Freeze()
        if not self.tags_root:
            self.tags_root = self.TreeTags.AddRoot(text="DICOM Objects")
        else:
            self.TreeTags.DeleteChildren(self.tags_root)
        self._recurse_tags_tree(ds, self.tags_root)
        self.TreeTags.ExpandAll()
        self.TreeTags.ScrollTo(self.TreeTags.GetFirstChild(self.tags_root)[0])
        self.Thaw()


    def _recurse_tags_tree(self, ds, parent, hide=False):
        """ order the dicom tags """
        for data_element in ds:
            if isinstance(data_element.value, pydicom.compat.text_type):
                text = pydicom.compat.text_type(data_element)
                ip = self.TreeTags.AppendItem(parent, text=text)
            else:
                ip = self.TreeTags.AppendItem(parent, text=str(data_element))

            if data_element.VR == "SQ":
                for i, ds in enumerate(data_element.value):
                    item_describe = data_element.name.replace(" Sequence", "")
                    item_text = "%s %d" % (item_describe, i + 1)
                    rjust = item_text.rjust(128)
                    parentNodeID = self.TreeTags.AppendItem(ip, text=rjust)
                    self._recurse_tags_tree(ds, parentNodeID)


    def _set_status(self, message):
        self.LabelCount.SetLabel(message)
        self.LabelCount.Update()

# ------------------------------------------------------------------------------
# Helper Classes for sorting data

class DicomFileInfo(object):
    """
    A convenience container for a small subset of the attributes in a DICOM file.

    The advantages of this over a pydicom Dataset object are --
    1) This uses Pythonic names
    2) All attributes are guaranteed to exist (although they may be empty)
       so there is no need to test "if 'foo' in dataset".
    3) This offers flexible sorting. Lists of DicomFileInfo objects will sort on
       the attributes listed in the class-level attribute 'sort_attributes'.
       Ex. To sort on 'study_id' and 'series_number':
           files = get_all_files(path)
           DicomFileInfo.sort_attributes = ("study_id", "series_number")
           files.sort()
           
    One cannot sort on arbitrary pydicom dataset objects, only on instance attributes.

    """
    sort_attributes = ("patient_name", "study_id", "series_number", "instance_number")

    def __init__(self, dataset, filename=""):
        """
        The dataset param must be a pydicom dataset. The filename is stored in
        this object for convenience but is otherwise unused.

        'slice_location' is calculated from ImagePositionPatient (0020 0037)
        and ImageOrientationPatient (0020 0037) if present according to the
        algorithm outlined at great length by Jolinda Smith here:
            http://www.itk.org/pipermail/insight-users/2003-September/004762.html
            http://www.cmake.org/pipermail/insight-users/2003-September/004762.html
        If those links die you can beseech Google:
            http://www.google.com/search?q=%22slightly-rotated+coronal+acquisition%22
        
        """
        self.dataset = dataset
        self.filename = filename
        
        self.patient_name       = str(dataset.PatientName)  # 0x0010 0010
        self.image_type         = 0                         # 0008 0008
        self.patient_id         = dataset.PatientID         # 0010 0020
        self.study_instance_uid = dataset.StudyInstanceUID  # 0020 000D
        self.study_id           = ""                        # 0020 0010
        self.study_description  = ""                        # 0008 1030
        self.image_id           = 0                         # 0054 0400
        self.instance_number    = 0                         # 0020 0013
        self.series_number      = 0                         # 0020 0011
        self.series_description = ""                        # 0008 103E
        self.series_instance_uid = ""                       # 0020 000E
        self.echo_time          = 0.001                     # 0018 0081
        self.acquisition_time   = 0                         # 0008 0032

        # Here I populate the not-always-present values
        if "SeriesInstanceUID" in dataset: self.series_instance_uid = dataset.SeriesInstanceUID
        if "SeriesDescription" in dataset: self.series_description  = dataset.SeriesDescription
        if "SeriesNumber"      in dataset: self.series_number       = dataset.SeriesNumber
        if "StudyID"           in dataset: self.study_id            = dataset.StudyID
        if "StudyDescription"  in dataset: self.study_description   = dataset.StudyDescription
        if "ImageID"           in dataset: self.image_id            = dataset.ImageID
        if "InstanceNumber"    in dataset: self.instance_number     = dataset.InstanceNumber
        if "EchoTime"          in dataset: self.echo_time           = dataset.EchoTime
        if "AcquisitionTime"   in dataset: self.acquisition_time    = dataset.AcquisitionTime
        if "ImageType"         in dataset: self.image_type          = dataset.ImageType

        self.slice_location = 0
            
        if  ("ImagePositionPatient" in dataset)     and \
            ("ImageOrientationPatient" in dataset)  and \
            dataset.ImagePositionPatient            and \
            (len(dataset.ImageOrientationPatient) >= 6):
           
            o = dataset.ImageOrientationPatient
    
            slice_normals = [ (o[1] * o[5]) - (o[2] * o[4]),
                              (o[2] * o[3]) - (o[0] * o[5]),
                              (o[0] * o[4]) - (o[1] * o[3]), ]
           
            self.slice_location = sum([a * b for a, b in zip(slice_normals, dataset.ImagePositionPatient)])
                            
    # bjs June 2020, Py3 no longer supports __cmp__() method, have to do all comparisons below

    def __lt__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 < item2

    def __gt__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 > item2

    def __eq__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 == item2

    def __ne__(self, other):
        """ Compares based on the class attribute sort_attributes"""
        item1, item2 = self.get_cmp_lists(other)
        return item1 != item2

    def get_cmp_lists(self, other):
        item1 = [getattr(self, item) for item in self.sort_attributes]
        item2 = [getattr(other, item) for item in self.sort_attributes]
        return item1, item2        

        
    def get_image(self):
        """
        Returns a wx.Image containing the image data in the dataset. 
        Only a small subset of the image types described in the DICOM 
        specification are supported. Specifically --
        - the dataset's PixelData attribute must be populated
        - the image must be grayscale
        - the data must be 8, 16, or 32-bit 
        - the transfer syntax UID must indicate raw mode (big and little
          endian are both OK)

        If this method encounters something it doesn't understand, it raises
        a NotImplementedError.
        """
        import wx
        
        ds = self.dataset       # shortcut ref
        
        if "PixelData" not in ds:
            return None
        
        else:
            RAW_UIDS = (pydicom.uid.ExplicitVRLittleEndian, 
                        pydicom.uid.ExplicitVRBigEndian)

            if ds.file_meta.TransferSyntaxUID not in RAW_UIDS:
                msg = "Unsupported: TransferSyntaxUID (0002,0010) == %d" % ds.file_meta.TransferSyntaxUID
                raise NotImplementedError(msg)

            if "PhotometricInterpretation" in ds:                       # We only handle monochrome images
                if not ds.PhotometricInterpretation.startswith("MONOCHROME"):
                    msg = "Unsupported: PhotometricInterpretation (0028,0004) == %d" % ds.PhotometricInterpretation
                    raise NotImplementedError(msg)

            data = ds.pixel_array
            nx, ny = data.shape

            if ("WindowCenter" in ds) and ("WindowWidth" in ds):
                cen = ds['WindowCenter']
                wid = ds['WindowWidth']
                cen = int(cen.value[0] if cen.VM > 1 else cen.value)
                wid = int(wid.value[0] if wid.VM > 1 else wid.value)
            else:
                wid = int(data.max() - data.min())
                cen = int(wid/2 + data.min())

            data = self.normalize_pixels(data, cen, wid)             # scale to 0-255

            if "PhotometricInterpretation" in ds:                   # MONOCHROME1 ==> lo vals bright / hi vals dark
                if ds.PhotometricInterpretation == "MONOCHROME1":
                    data = 255 - data

            # wx.Image wants RGB byte array
            data = data.astype(np.uint8)        # may not be necessary, bjs ?
            data = np.repeat(data, 3)
            data = data.tobytes()

            image = wx.Image((nx, ny), data)

            return image


    def normalize_pixels(self, data, _window_center, _window_width):
        """
        Applies a linear conversion to the pixels based on the Window
        Center (0028,1050) and Window Width (0028,1051) values of the dataset.
        Data is a numpy array. It is altered in-place and returns floats
        between 0 and 255.

        This implements the algorithm shown in the DICOM specification
        Part 3: Information Object Definitions, Annex C.11.2.1.2.
        ftp://medical.nema.org/medical/dicom/2008/08_03pu.pdf

        Here's the pseudocode from the DICOM spec:

            if      (x <= c - 0.5 - (w-1)/2), then y = ymin
            else if (x >  c - 0.5 + (w-1)/2), then y = ymax
            else y = ((x - (c - 0.5)) / (w-1) + 0.5) * (ymax - ymin) + ymin

        """
        win = int(_window_width)
        lvl = int(_window_center)
        e = [0, 255, lambda data: ((data - (lvl - 0.5)) / (win - 1) + 0.5) * (255 - 0)]
        return np.piecewise(data, [data <= (lvl - 0.5 - (win - 1) / 2), data > (lvl - 0.5 + (win - 1) / 2)], e)


class DicomFileSiemens(DicomFileInfo):
    """
    A DICOM file container with some Siemens-specific qualities. This
    is not anywhere close to a full implementation of Siemens' DICOM 
    extensions. Rather, it handles a couple of specific things we care about.
    
    """
    TAG_CONTENT_TYPE = (0x0029, 0x1008)
    TAG_IMAGE_TYPE = (0x0008, 0x0008)

    # examples
    # CS: ['DERIVED', 'SECONDARY', 'MPR', 'CSA MPR', '', 'CSAPARALLEL', 'M', 'ND']
    # CS: ['ORIGINAL', 'PRIMARY']
    # CS: ['DERIVED', 'SECONDARY', 'SPEC']

    def __init__(self, dataset, filename=""):
        
        DicomFileInfo.__init__(self, dataset, filename)
        
        self.content_type = TYPE_NONE
        content_type = ""
        
        # pydicom doesn't know about the Siemens content type tag, so I have 
        # to pass the magic numbers that represent it.
        if self.TAG_CONTENT_TYPE in dataset:
            content_type = dataset[self.TAG_CONTENT_TYPE].value.upper()
            if "SPEC" in content_type:
                self.content_type = TYPE_SPECTROSCOPY
            if "IMAGE" in content_type:
                self.content_type = TYPE_IMAGE
        else:
            if self.TAG_IMAGE_TYPE in dataset:
                content_type = dataset[self.TAG_IMAGE_TYPE].value
                if "MPR" in content_type and "DERIVED" in content_type:
                    self.content_type = TYPE_IMAGE
            else:
                content_type = ""
        
        # There are lots of values that can appear in the Siemens content 
        # type tag. The two we care about start with SPEC and IMAGE.
        # ref: http://www.medical.siemens.com/webapp/wcs/stores/servlet/CategoryDisplay~q_catalogId~e_-11~a_categoryId~e_16560~a_catTree~e_100003,16554,16560~a_langId~e_-11~a_storeId~e_10001.htm
        # alt ref: http://www.medical.siemens.com/DICOM

        self.dataset = dataset



# =============================================================================
# Dialog DICOM Browser Subclasses that show Manufacturer Specific Files

class SiemensMrsBrowser(DialogDicomBrowser):
    """
    This sub-class will select for Siemens and Spectroscopy files ONLY

    """
    TAG_SOP_CLASS_UID = (0x0008, 0x0016)
    TAG_CONTENT_TYPE = (0x0029, 0x1008)
    TAG_IMAGE_TYPE = (0x0008, 0x0008)


    def __init__(self, parent=None,
                       preview_size=(64, 64),
                       show_tags=None,
                       multi_select=False,
                       return_objects=False,
                       default_path=''):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        self.manufacturer = ['SIEMENS', ]
        DialogDicomBrowser.__init__(self, parent,
                                          preview_size=preview_size,
                                          show_tags=show_tags,
                                          multi_select=multi_select,
                                          return_objects=return_objects,
                                          default_path=default_path)

    def file_filter(self, dicom_file):
        """
        This function is called for each file to see whether or not it should
        be included in the tree. It returns True if the given DICOM file is a
        a Siemens spectroscopy file.

        """
        dataset = dicom_file.dataset
        test = False
        if self.manufacturer:
            if 'Manufacturer' in dataset:
                for item in self.manufacturer:
                    # we only handle manufacturers in this list
                    if dataset.Manufacturer.upper() in item:
                        test = True
            if not test:
                return False

        if self.TAG_SOP_CLASS_UID in dataset:
            item_uid = dataset[self.TAG_SOP_CLASS_UID].value.upper()
            sop_class_uid = pydicom.uid.UID(str(item_uid))
            if sop_class_uid.name == 'MR Image Storage':
                # this is the condition for Spectroscopy Browser
                return False
            elif sop_class_uid.name == 'MR Spectroscopy Storage':
                # this is the UID value for Siemens VD files
                return True
            elif sop_class_uid.name == str(item_uid):
                # In software versions VA, VB (maybe others, but not VD), Siemens uses a proprietary
                # SOP Class UID for their spectroscopy data files. I have seen 1.3.12.2.1107.5.9.1
                # used but am not sure if it is unique. And this is not in the DICOM dictionary of UIDs
                # So, we will look for other proprietary tags that will give us solid info that these are
                # Siemens SPEC file.
                #
                # pydicom doesn't know about the Siemens content type tag, so I have
                # to pass the magic numbers that represent it.
                if self.TAG_CONTENT_TYPE in dataset:
                    content_type = dataset[self.TAG_CONTENT_TYPE].value.upper()
                else:
                    content_type = ""
                # There are lots of values that can appear in the Siemens content
                # type tag. The ones I care about are "SPEC NUM 4" and "Spectroscopy".
                # ref:
                #    p 132 of this PDF:
                #    http://www.medical.siemens.com/siemens/en_GLOBAL/rg_marcom_FBAs/files/brochures/DICOM/mr/syngo_MR_B17.pdf
                #    which is from here:
                #    http://www.medical.siemens.com/webapp/wcs/stores/servlet/CategoryDisplay~q_catalogId~e_-11~a_categoryId~e_16560~a_catTree~e_100003,16554,16560~a_langId~e_-11~a_storeId~e_10001.htm

                return ("SPEC" in content_type)
            else:
                return False

class PhilipsMrsBrowser(DialogDicomBrowser):
    """
    This sub-class will select for Philips and Spectroscopy files ONLY

    """

    def __init__(self, parent=None,
                       preview_size=(64, 64),
                       show_tags=None,
                       multi_select=False,
                       return_objects=False,
                       default_path=''):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        self.manufacturer = ['PHILIPS', ]
        DialogDicomBrowser.__init__(self, parent,
                                          preview_size=preview_size,
                                          show_tags=show_tags,
                                          multi_select=multi_select,
                                          return_objects=return_objects,
                                          default_path=default_path)

    def file_filter(self, dicom_file):
        """
        This function is called for each file to see whether or not it should
        be included in the tree. It returns True if the given DICOM file is a
        a Philips spectroscopy file.

        """
        dataset = dicom_file.dataset
        test = False
        if self.manufacturer:
            if 'Manufacturer' in dataset:
                for item in self.manufacturer:
                    # we only handle manufacturers in this list
                    if dataset.Manufacturer.upper() in item:
                        test = True
            if not test:
                return False

        return self.is_mrs_dicom(dicom_file.dataset)


    def is_mrs_dicom(self, dataset):
        """ returns True if all criteria (Sandeep) are met """

        if not isinstance(dataset, pydicom.dataset.Dataset):
            raise ValueError("Object passed in not a dicom Dataset.")

        if not "ProtocolName" in dataset:
            return False

        if dataset.ProtocolName == 'ExamCard':
            return False

        if not (0x2005, 0x10c0) in dataset:
            return False

        if dataset[0x2005, 0x10c0].value != 'SPECTRO':  # TAG_PHILIPS_PRIVATE_02
            return False

        if not (0x2005, 0x1270) in dataset:             # TAG_PHILIPS_PRIVATE_01
            return False

        return True


# -----------------------------------------------------------------------------
# Embedded Icons

book_icon = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABGdBTUEAAK/INwWK6QAAABl0'
    b'RVh0U29mdHdhcmUAQWRvYmUgSW1hZ2VSZWFkeXHJZTwAAAHjSURBVDjLdZO/alVBEMZ/5+Te'
    b'mxAbFUUskqAoSOJNp4KC4AsoPoGFIHY+gA+jiJXaKIiChbETtBYLUbSMRf6Aydndmfks9kRj'
    b'vHdhGVh2fvN9uzONJK7fe7Ai6algA3FZCAmQqEF/dnihpK1v7x7dPw0woF64Izg3Xl5s1n9u'
    b'Ie0lQYUFCtjc+sVuEqHBKfpVAXB1vLzQXFtdYPHkGFUCoahVo1Y/fnie+bkBV27c5R8A0pHx'
    b'yhKvPn5hY2MHRQAQeyokFGJze4cuZfav3gLNYDTg7Pklzpw4ijtIQYRwFx6BhdjtCk+erU0C'
    b'CPfg+/o2o3ZI13WUlLGo58YMg+GIY4dmCWkCAAgPzAspJW5ePFPlV3VI4uHbz5S5IQfy/yoo'
    b'HngxzFser30iFcNcuAVGw3A0Ilt91IkAsyCXQg5QO0szHEIrogkiguwN2acCoJhjnZGKYx4U'
    b'jz5WOA2YD1BMU+BBSYVUvNpxkXuIuWgbsOxTHrG3UHIFWIhsgXtQQpTizNBS5jXZQkhkcywZ'
    b'qQQlAjdRwiml7wU5xWLaL1AvZa8WIjALzIRZ7YVWDW5CiIj48Z8F2pYLl1ZR0+AuzEX0UX03'
    b'5mxIkLq0dhDw5vXL97fr5O3rfwQHJhPx4uuH57f2AL8BfPrVlrs6xwsAAAAASUVORK5CYIIj'
    b'LS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0t'
    b'LS0tLS0tLS0tLS0tLS0tLQ0KZmxhdGljb25fMzYzOF8xNiA9IFB5RW1iZWRkZWRJbWFnZSgN'
    b'CiAgICBiJ2lWQk9SdzBLR2dvQUFBQU5TVWhFVWdBQUFCQUFBQUFRQ0FZQUFBQWY4LzloQUFB'
    b'QUJITkNTVlFJQ0FnSWZBaGtpQUFBQUFsdycNCiAgICBiJ1NGbHpBQUFBZGdBQUFIWUJUbnNt'
    b'Q0FBQUFCbDBSVmgwVTI5bWRIZGhjbVVBZDNkM0xtbHVhM05qWVhCbExtOXlaNXZ1UEJvQScN'
    b'CiAgICBiJ0FBRVZTVVJCVkRpTm5kSy9Ma1JSRUFid0g1Rm8yR3h0TldvcU5CSUZpWVpFTC9F'
    b'U1BJTkN5VU1RSG1MOXJVaFdwNVFOaVYydCcNCiAgICBiJ0NJWEdLbmJ1N3RtOUo0UXZtWHoz'
    b'emprelorYWJvWXdLN3REQmZsZ0hqVGdid0dnbXdRYm1NLzRGYkE3SHBEODEzT01rRTF6ZycN'
    b'CiAgICBiJ0NNOVlMQnhqeWVFT1pxUFUxL0E5Qk5lRHF4RzhoL1hoQkZQQngxRUpQQVkvQmM5'
    b'Rmd0cHdhWk00MXhVcnRWVEUxSnFZb0s5QicNCiAgICBiJ0hhdjV0ck9ZS2RvYWlkZmY0dUFR'
    b'bjhuRnkrQ1Z4RGV1cXhkVWNtUDhGMjZVKy94SmcwN0U5RFJZdzhVZkhteEdURytNNzNpSicN'
    b'CiAgICBiJzcxM2xNWjRGeitFQUh4RXpzQWZ0NEczOVJUb04zZ3F1QnJkeVpSV3IvSnNHYmNr'
    b'cXAxTm9SWWxMdWV5QlpVenJybnNwUVlFRycNCiAgICBiJ3JqUCtLOXppSzNWK0E1N0hUSFo1'
    b'RVJMVUFBQUFBRWxGVGtTdVFtQ0MnKQ0KDQo=')

#----------------------------------------------------------------------
top_icon = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABHNCSVQICAgIfAhkiAAAAAlw'
    b'SFlzAAAAdgAAAHYBTnsmCAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoA'
    b'AAEVSURBVDiNndK/LkRREAbwH5Fo2GxtNWoqNBIFiYZEL/ESPINCyUMQHmL9rUhWp5QNiV2t'
    b'CIXGKnbu7tm9J4QvmXz3zjkzZ+aboYwK7tDBflgHjTgbwGgmwQbmM/4FbA7HpD813OMkE1zg'
    b'CM9YLBxjyeEOZqPU1/A9BNeDqxG8h/XhBFPBx1EJPAY/Bc9FgtpwaZM41xUrtVTE1JqYoK9B'
    b'Hav5trOYKdoaidff4uAQn8nFy+CVxDeuqxdUcmP8F26U+/xJg07E9DRYw8UfHmxGTG+M73iJ'
    b'713lMZ4Fz+EAHxEzsAft4G39RToN3gquBrdyZRWr/JsGbckqp1NoRYlLueyBZUzrrnspQYEG'
    b'rjP+K9ziK3V+A57HTHZ5ERLUAAAAAElFTkSuQmCC')

#----------------------------------------------------------------------
mri_icon = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAAAAXNSR0IArs4c6QAAAARnQU1B'
    b'AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAZdEVYdFNvZnR3YXJlAHBhaW50Lm5l'
    b'dCA0LjAuMTnU1rJkAAAArUlEQVQ4T2P4jw00d/bv3ncIykEFKBps3QIwkZ17AFQaDBAa0NSh'
    b'IagiuAY0aawIohKkAVMUDiITs+BSh4+dBIqgaACrQQdwWYgCFE/jAigaLGsOQhCQg8xGBtg1'
    b'ICOIOjjArgEih6kaCFA0eLUfw4ogSoHA0SsERQNEFBeYPnsBXPXCpauAImRFHASgSaMhqCI0'
    b'J6EpgiOoNBhg8cO/f/88AqIy8suhfGTw/z8AaPkEJ2t0SU0AAAAASUVORK5CYII=')

#----------------------------------------------------------------------
user_icon = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABGdBTUEAAK/INwWK6QAAABl0'
    b'RVh0U29mdHdhcmUAQWRvYmUgSW1hZ2VSZWFkeXHJZTwAAAJ3SURBVDjLpZNtSNNRFIcNKunF'
    b'1rZWBMJqKaSiX9RP1dClsjldA42slW0q5oxZiuHrlqllLayoaJa2jbm1Lc3QUZpKFmmaTMsa'
    b'Rp+kMgjBheSmTL2//kqMBJlFHx44XM7vOfdyuH4A/P6HFQ9zo7cpa/mM6RvCrVDzaVDy6C5J'
    b'JKv6rwSnIhlFd0R0Up/GwF2KWyl01CTSkM/dQoQRzAurCjRCGnRUUE2FaoSL0HExiYVzsQwc'
    b'j6RNrSqo4W5Gh6Yc4+1qDDTkIy+GhYK4nTgdz0H2PrrHUJzs71NQn86enPn+CVN9GnzruoYR'
    b'63mMPbkC59gQzDl7pt7rc9f7FNyUhPY6Bx9gwt4E9zszhWWpdg6ZcS8j3O7zCTuEpnXB+3MN'
    b'ZkUUZu0NmHE8XsL91oSWwiiEc3MeseLrN6woYCWa/Zl8ozyQ3w3Hl2lYy0SwlCUvsVi/Gv2J'
    b'wITnYPDun2Hy6jYuEzAF1jUBCVYpO6kXo+NuGMeBAgcgfwNkvgBOPgUqXgKvP7rBFvRhE1cr'
    b'p8Vq1noFYSlacVyqGk0D86gbART9BDk9BFnPCNJbCY5aCFL1Cyhtp0RWAp74MsKSrkq9guHy'
    b'vfMTtmLc1togpZoyqYmyNoITzVTYRJCiXYBIQ3CwFqi83o3JDhX6C0M8XsGIMoQ4OyuRlq1D'
    b'dZcLkmbgGDX1iIEKNxAcbgTEOqC4ZRaJ6Ub86K7CYFEo8Qo+GBQlQyXBczLZpbloaQ9k1NUz'
    b'/kD2myBBKxRZpa5hVcQslalatoUxizxAVVrN3CW21bFj9F858Q9dnIRmDyeuybM71uxmH9BN'
    b'BB1q6zybV7H9s1Ue4PM3/gu/AEbfqfWy2twsAAAAAElFTkSuQmCC')






# ------------------------------------------------------------------------------
# Test routines


class TestFrame(wx.Frame):
    """This frame demonstrates calling the DICOM file browser dialog."""

    msg = """
This is a dummy window to provide a parent for the DICOM file
browser dialog. In a real application, this window would be 
replaced by the application's main window.
"""

    def __init__(self, default_path=''):
        wx.Frame.__init__(self, None, title="DICOM Viewer", size=(500, 130))

        label = wx.StaticText(self, -1, self.msg)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(label, 0, wx.ALL, 10)
        self.SetSizer(sizer)
        self.SetBackgroundColour("white")   # in Windows, background gray if unspecified
        self.Layout()
        self.Center()
        self.Show()

        # preview = None (default), disable preview feature
        # show_tags = None (default), disable checkbox 'show tags' feature

        # dialog = DialogDicomBrowser(self,
        #                             preview_size=(96, 96),
        #                             show_tags=False,
        #                             multi_select=True,
        #                             return_objects=False,
        #                             default_path=default_path)

        # dialog = SiemensSpectroscopyDicomBrowserDialog(self,
        #                             preview_size=(96, 96),
        #                             show_tags=False,
        #                             multi_select=True,
        #                             return_objects=False,
        #                             default_path=default_path)

        dialog = PhilipsSpectroscopyDicomBrowserDialog(self,
                                    preview_size=(96, 96),
                                    show_tags=False,
                                    multi_select=True,
                                    return_objects=False,
                                    default_path=default_path)
        dialog.ShowModal()

        if dialog.items:
            print("You selected %d item(s) --" % len(dialog.items))
            for item in dialog.items:
                if hasattr(item, "filename"):
                    print('- ' + item.filename)
                else:
                    print('- ' + item)
        else:
            print("You opted not to select an item.")

        print("\nbye!")

        self.Close()


def _test():
    default_path = ''
    default_path = r'D:\Users\bsoher\code\repository_svn\sample_data\siemens_dicom_export'

    app = wx.App(0)
    wx.InitAllImageHandlers()
    frame_1 = TestFrame(default_path=default_path)
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()


if __name__ == '__main__':
    _test()
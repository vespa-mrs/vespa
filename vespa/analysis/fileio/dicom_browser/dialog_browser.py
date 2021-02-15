# Python modules

import os.path

# 3rd party modules
import wx
try:
    import dicom
except:
    import pydicom as dicom

# Local modules
import vespa.analysis.fileio.dicom_browser.auto_gui.browser as browser
import vespa.common.wx_gravy.util as wx_util
import vespa.common.util.misc as util_misc
import vespa.common.util.config as util_config
import vespa.analysis.fileio.dicom_browser.util_browser as util_browser
import vespa.analysis.fileio.util_philips as util_philips

TAG_SOP_CLASS_UID = (0x0008, 0x0016)
# This tag holds a description of the content of a Siemens file
TAG_CONTENT_TYPE  = (0x0029, 0x1008)


class BrowserDialog(browser.TheDialog):
    """Opens a dialog that allows the user to select a directory containing
    DICOM files. Displays the files in a tree and allows the user to select
    one.

    By default a preview is displayed if the DICOM file contains an image
    that our code understands.

    After the dialog closes, the selected filenames are stored in the
    filenames attribute of this class. Note that for consistency's sake, 
    dialog.filenames is a list even if the dialog is invoked with the 
    multi_select option set to False.
    
    If the user doesn't select a file, dialog.filenames is an empty list.
    """

    def __init__(self, parent=None, preview_size=(64, 64), multi_select=False):
        """Creates an instance of the dialog.

        If preview_size is None, the dialog won't display previews of
        image files. Otherwise the preview size should be a tuple that
        describes the preview size.

        If multi_select is False, the tree only allows single selection.
        """
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        style = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        browser.TheDialog.__init__(self, parent, -1, style=style)

        self.multi_select = multi_select

        if multi_select:
            style = self.Tree.GetWindowStyle()
            style &= ~wx.TR_SINGLE
            style |= wx.TR_MULTIPLE 
            #tyle |= (wx.TR_MULTIPLE | wx.TR_EXTENDED)
            self.Tree.SetWindowStyle(style)
            
            
        # We add the Open & Cancel buttons dynamically so that they're in  
        # the right order under OS X, GTK, Windows, etc.
        self.ButtonOpen, _ = wx_util.add_ok_cancel(self, 
                                                   self.LabelOpenCancelPlaceholder, 
                                                   self.on_open, 
                                                   ok_text="Open")
        
        if preview_size:
            # Add the preview bitmap as the leftmost item in the grid sizer.
            self.preview = wx.StaticBitmap(self, -1, wx.NullBitmap)
            self.preview.SetMinSize(preview_size)
            grid_sizer.Insert(0, self.preview, 0, wx.LEFT|wx.BOTTOM, 10)
        else:
            # Caller doesn't want a preview
            self.preview = None

        self.SetSize( (700, 350) )

        self.Layout()

        # I add a fake root to the tree. The fake root will be invisible. This
        # allows me to display multiple sets of DICOM data as siblings.
        self.Tree.AddRoot("")

        # These must be set before anything reads from or writes to the
        # property self.path.
        self._path = ""
        self._user_home = wx.GetUserHome()

        self.path = ""

        self.filenames = [ ]
        
        self._set_status("")

        self._eat_a_selection_event = False

        # Here I set up the icon used in the tree.
        image_size = (16, 16)
        images = wx.ImageList(image_size[0], image_size[1])

        the_icon = wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, image_size)
        self.icon_index = images.Add(the_icon)
        self.Tree.AssignImageList(images)

        self.Center()

        self.Bind(wx.EVT_TREE_SEL_CHANGING, self.on_tree_selection_changing, self.Tree)

    __doc = """User's currently selected path which is reflected in the
    dialog's title bar."""
    def __get_path(self):
        return self._path
    def __set_path(self, path):
        self._path = path
        if not path:
            path = "[no directory selected]"
        elif path.startswith(self._user_home):
            # Trim off the home dir + separator to save space
            path = path[len(self._user_home) + 1:]

        self.SetTitle("DICOM Browser - " + path)
    path = property(__get_path, __set_path, doc=__doc)


    ####################    Public methods    #####################
    
    def get_item_description(self, dicom_file):
        """Builds the description string for a leaf item in the list and
        returns that string along with an index pointing to the appropriate
        icon for this item in the tree's ImageList.

        If you want a different description, subclass this dialog and
        override this function. Return -1 for the icon index if you don't
        want an icon.
        """
        description = "%d %s" % (dicom_file.instance_number,
                                 os.path.split(dicom_file.filename)[1])

        return description, self.icon_index

    
    def file_filter(self, dicom_file):
        """This function is called for each file to see whether or not it 
        should be included in the tree. If the function returns False, the
        file is filtered out (excluded).

        This version of the function excludes no files because it always
        returns True. It exists so that one can implement filtering by 
        subclassing the dialog and overriding this function.
        """
        return True


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
            if self.Tree.ItemHasChildren(item_id):
                # Select all kids. Selecting items tends to scroll the list. We 
                # save the first item so we can restore the current scroll 
                # position when selecting is done. We also disable selection
                # notifications while doing this -- it speeds things up 
                # considerably.
                first_visible = self.Tree.GetFirstVisibleItem()
                self.Unbind(wx.EVT_TREE_SEL_CHANGED, self.Tree)
                self._select_offspring(item_id)
                self.Bind(wx.EVT_TREE_SEL_CHANGED, self.on_tree_selection_changed, self.Tree)

                # Strange but true: under OS X & GTK, an EVT_LEFT_UP event 
                # occurs after this event, and that triggers an item selection 
                # which undoes the work of _select_offspring(). To work around 
                # this, I set a flag which eats the next selection event. 
                if wx.Platform != "__WXMSW__":
                    self._eat_a_selection_event = True

                self.Tree.ScrollTo(first_visible)
        # else:
        #     If only single-selection is enabled, or multiple items have been
        #     selected, don't enable the header-implies-selection behavior.


        filenames = self._get_selected_filenames()
        
        # The Open button is enabled as long as at least one filename
        # is selected.
        self.ButtonOpen.Enable(bool(filenames))

        self._set_status("{0} of {1} files selected.".format(len(filenames),
                                                             self._nfiles))

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
                item = self.Tree.GetItemData(item_ids[0])
                if item:
                    filename = item.filename

                    try:
                        image = item.get_image()
                    except NotImplementedError:
                        # This DICOM file contains an image in a format our
                        # code doesn't understand.
                        image = None

                    if image:
                        size = self.preview.GetSize()
                        image.Rescale(size[0], size[1])
                        bitmap = wx.BitmapFromImage(image)
                #else:
                    # The item selected is not a leaf node
            #else:
                # There's more than one item selected

            self.preview.SetBitmap(bitmap)


    def on_button_choose_directory(self, event):
        self.path = wx.DirSelector("Select a directory of DICOM files",
                                   self._suggest_a_start_path())

        if self.path:
            self._populate_tree()


    def on_tree_item_activated(self, event):
        # "Activated" == a double click
        # I don't want a double click to do anything if the currently
        # selected item isn't a leaf node.
        item_id = event.GetItem()

        if item_id and (not self.Tree.ItemHasChildren(item_id)):
            self.on_open(event)


    def on_open(self, event):
        self.filenames = self._get_selected_filenames()

        # Record the current path
        if self.filenames:
            path, _ = os.path.split(self.filenames[0])

            config = util_config.VespaConfig()
            config["general"]["last_dicom_browse_path"] = path
            config.write()

        self.Close()


    ####################    "Private" methods    #####################

    def _get_selections(self):
        if self.multi_select:
            item_ids = self.Tree.GetSelections()
        else:
            item_ids = [ self.Tree.GetSelection() ]

        return item_ids


    def _get_selected_filenames(self):
        item_ids = self._get_selections()

        items = [self.Tree.GetItemData(item_id) for item_id in item_ids]
        
        # The filenames need to be converted to Unicode. See comment in 
        # _populate_tree() regarding pydicom issue 58.
        return [str(item.filename) for item in items if item]


    def _select_offspring(self, item_id):
        """Given a tree item id, recursively selects that item's children, 
        including grandchildren, great-grandchildren, etc.
        """
        kid_id, cookie = self.Tree.GetFirstChild(item_id)
        while kid_id.IsOk():
            self.Tree.SelectItem(kid_id)
            if self.Tree.ItemHasChildren(kid_id):
                self._select_offspring(kid_id)
            kid_id, cookie = self.Tree.GetNextChild(item_id, cookie)


    def _set_status(self, message):
        self.LabelCount.SetLabel(message)
        self.LabelCount.Update()


    def _suggest_a_start_path(self):
        """Returns a path for the user to begin browsing for DICOM files."""
        # If the user already selected a path, I start there.
        path = self.path

        if not os.path.exists(path):
            # Maybe there's one saved in the config file
            config = util_config.VespaConfig()
            section = config["general"]
            path = section.get("last_dicom_browse_path", "")   

        if not os.path.exists(path):
            path = util_misc.get_documents_dir()

        return path


    def _populate_tree(self):
        wx.SetCursor(wx.HOURGLASS_CURSOR)

        root = self.Tree.GetRootItem()

        # Reset controls
        self.Tree.DeleteChildren(root)
        self.ButtonOpen.Enable(False)
        self._set_status("Scanning for DICOM files...")
        self.filenames = [ ]
        if self.preview:
            self.preview.SetBitmap(wx.NullBitmap)

        # Under Windows & GTK, this dialog needs a repaint after the file
        # open dialog disappears. Unless I force one manually here, it won't
        # happen until this lengthy function completes. It doesn't always
        # get the job done under GTK. :-/
        self.Update()

        # bjs June 2020, not needed in Py3
        # # Per issue 58, pydicom doesn't handle paths that are Unicode strings
        # # properly (even if the contents are 100% ASCII), so I make sure the
        # # path I pass is a string.
        # # http://code.google.com/p/pydicom/issues/detail?id=58&can=1&q=defer_size&colspec=ID%20Type%20Status%20Priority%20Milestone%20Owner%20Summary%20Difficulty
        # path = self.path
        # if isinstance(self.path, str):
        #     path = path.encode("utf-8")

        files = [ ]
        self._nfiles = 0
        for f in util_browser.get_files(self.path):
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
        root = self.Tree.GetRootItem()
        for df in files:
            if df.patient_id != patient_id:
                # New patient
                study_id = None
                series_number = None
                patient_id = df.patient_id
                patient_name = df.patient_name
                s = "%s : %s" % (df.patient_id, df.patient_name)
                current_patient = self.Tree.AppendItem(root, s)

            if df.study_id != study_id:
                # New study
                series_number = None
                study_id = df.study_id
                study_description = df.study_description
                s = "Study %s %s" % (df.study_id, df.study_description)
                current_study = self.Tree.AppendItem(current_patient, s)

            if df.series_number != series_number:
                series_number = df.series_number
                series_description = df.series_description
                s = "Series %d %s" % (df.series_number, df.series_description)
                current_series = self.Tree.AppendItem(current_study, s)

            description, icon_index = self.get_item_description(df)

            item_id = self.Tree.AppendItem(current_series, description,
                                           icon_index)
            # Each leaf node's data is set to the associated DicomFileInfo instance.
            self.Tree.SetItemData(item_id, df)

        self.Tree.ExpandAll()

        if files:
            # Under Windows the list is scrolled to the bottom, I think 
            # because of the call to ExpandAll() above. We force the first
            # item to be visible.
            first_child, _ = self.Tree.GetFirstChild(root)
            self.Tree.ScrollTo(first_child)
        else:
            self._set_status("No DICOM files found.")

        wx.SetCursor(wx.NullCursor)


class SiemensSpectroscopyBrowserDialog(BrowserDialog):
    def __init__(self, parent=None, multi_select=False):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        BrowserDialog.__init__(self, parent, None, multi_select)
        
        self.manufacturer = ['SIEMENS',]

    def file_filter(self, dicom_file):
        """
        This function is called for each file to see whether or not it should
        be included in the tree. It returns True if the given DICOM file is a
        a Siemens spectroscopy file.
        
        """
        test = False
        
        dataset = dicom_file.dataset
        
        if self.manufacturer:
            if 'Manufacturer' in dataset:
                for item in manufacturer:
                    # we only handle manufacturers in this list
                    if dataset.Manufacturer.upper() in item:
                        test = True
            if not test:
                return False
        
        if TAG_SOP_CLASS_UID in dataset:
            item_uid = dataset[TAG_SOP_CLASS_UID].value.upper() 
            sop_class_uid = dicom.uid.UID(str(item_uid))
            if sop_class_uid.name == 'MR Image Storage':
                # this is the condition for Spectroscopy Browser
                return False
            elif sop_class_uid.name == 'MR Spectroscopy Storage':
                # this is the UID value for Siemens VD files
                return True
            elif sop_class_uid.name == str(item_uid):
                # In software versions VA, VB (maybe others, but not VD),
                # Siemens uses a proprietary SOP Class UID for their spectroscopy
                # data files. I have seen 1.3.12.2.1107.5.9.1 used but am not sure
                # if it is unique. And this is not in the DICOM dictionary of UIDs
                # So, we will look for other proprietary tags that will give us 
                # solid info that these are Siemens SPEC file.
                # 
                # pydicom doesn't know about the Siemens content type tag, so I have 
                # to pass the magic numbers that represent it.
                if TAG_CONTENT_TYPE in dataset:
                    content_type = dataset[TAG_CONTENT_TYPE].value.upper()
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
    


class PhilipsSpectroscopyBrowserDialog(BrowserDialog):
    def __init__(self, parent=None, multi_select=False):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        BrowserDialog.__init__(self, parent, None, multi_select)
        
        self.manufacturer=['PHILIPS',]


    def file_filter(self, dicom_file):
        """
        This function is called for each file to see whether or not it should
        be included in the tree. It returns True if the given DICOM file is a
        a Philips spectroscopy file.
        
        """
        test = False
        
        dataset = dicom_file.dataset

        if self.manufacturer:
            if 'Manufacturer' in dataset:
                for item in manufacturer:
                    # we only handle manufacturers in this list
                    if dataset.Manufacturer.upper() in item:
                        test = True
            if not test:
                return False
        
        return util_philips.is_mrs_dicom(dataset)  

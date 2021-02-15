# Python modules
import os
import os.path

# 3rd party modules
import wx
import pydicom 

# Our modules
import vespa.analysis.fileio.dicom_browser.auto_gui.dicom_browser as dicom_browser
import vespa.common.wx_gravy.util as wx_util
import vespa.analysis.fileio.util_philips as util_philips


TAG_SOP_CLASS_UID = (0x0008, 0x0016)
# This tag holds a description of the content of a Siemens file
TAG_CONTENT_TYPE  = (0x0029, 0x1008)


# MIN & MAX PIXEL are limited by what can be stuffed into a byte, because
# we turn the DICOM data into RGB data with R == G == B for each pixel.
_MAX_PIXEL = 255
_MIN_PIXEL = 0



class DicomBrowserDialog(dicom_browser.TheDialog):
    """
    Opens a dialog that allows the user to select a directory containing
    DICOM files. Displays files in a tree and allows the user to select one.

    By default a preview is displayed if the DICOM file contains an image
    that our code understands.

    After the dialog closes, the selected filenames are stored in the
    filenames attribute of this class. Note that for consistency's sake, 
    dialog.filenames is a list even if the dialog is invoked with the 
    multi_select option set to False.
    
    If the user doesn't select a file, dialog.filenames is an empty list.
    """

    def __init__(self, parent=None,
                       preview_size=(64, 64),
                       multi_select=False,
                       path=''):
        """Creates an instance of the dialog.

        If preview_size is None, the dialog won't display previews of
        image files. Otherwise the preview size should be a tuple that
        describes the preview size.

        If multi_select is False, the tree only allows single selection.
        """
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        style = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        dicom_browser.TheDialog.__init__(self, parent, -1, style=style)

        self.multi_select = multi_select

        if multi_select:
            style = self.Tree.GetWindowStyle()
            style &= ~wx.TR_SINGLE
            style |= wx.TR_MULTIPLE 
            self.Tree.SetWindowStyle(style)
            
            
        # dynam add Open & Cancel for right order under OSX, GTK, Windows, etc.
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
        self.Tree.AddRoot('')

        # Set these before anything reads/writes to the property self.path.
        self._path = path
        self._user_home = wx.GetUserHome()

        self.path = path
        self.filenames = []
        
        self._set_status('')

        self._eat_a_selection_event = False     # req for a bugfix, see below

        # Here I set up the icon used in the tree.
        image_size = (16, 16)
        images = wx.ImageList(image_size[0], image_size[1])
        the_icon = wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, image_size)
        self.icon_index = images.Add(the_icon)
        self.Tree.AssignImageList(images)

        self.Center()

        self.Bind(wx.EVT_TREE_SEL_CHANGING, self.on_tree_selection_changing, self.Tree)

        if self.path:
            self._populate_tree()
            self.TextCurrentDirectory.SetValue(self.path)

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


    def on_button_browse(self, event):
        """ Select a directory that has DICOM files """
        self.path = wx.DirSelector("Select a directory of DICOM files", self.path)
        if self.path:
            self.TextCurrentDirectory.SetValue(self.path)
            self._populate_tree()


    def on_tree_item_activated(self, event):
        """ 'Activated' == a double click only on a leaf node. """
        item_id = event.GetItem()

        if item_id and (not self.Tree.ItemHasChildren(item_id)):
            self.on_open(event)


    def on_open(self, event):
        """ Get filenames of items selected in tree """
        self.filenames = self._get_selected_filenames()
        if self.filenames:
            self.path, _ = os.path.split(self.filenames[0])  # store the current path
            #self.path = path
        self.Close()


    ####################    "Private" methods    #####################

    def _get_selections(self):
        """ return one or more selected wx object ids from tree """
        item_ids = self.Tree.GetSelections() if self.multi_select else [self.Tree.GetSelection(),]
        return item_ids


    def _get_selected_filenames(self):
        """ parse filenames from wx objects selected in the tree """
        items = [self.Tree.GetItemData(item_id) for item_id in self._get_selections()]
        return [str(item.filename) for item in items if item]


    def _select_offspring(self, item_id):
        """Given tree item id, recursively select all child, grandchild, etc items """
        kid_id, cookie = self.Tree.GetFirstChild(item_id)
        while kid_id.IsOk():
            self.Tree.SelectItem(kid_id)
            if self.Tree.ItemHasChildren(kid_id):
                self._select_offspring(kid_id)
            kid_id, cookie = self.Tree.GetNextChild(item_id, cookie)


    def _set_status(self, message):
        self.LabelCount.SetLabel(message)
        self.LabelCount.Update()


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

        files = [ ]
        self._nfiles = 0
        for f in get_files(self.path):
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
            if df.patient_id != patient_id:         # New patient
                study_id = None
                series_number = None
                patient_id = df.patient_id
                patient_name = df.patient_name
                s = "%s : %s" % (df.patient_id, df.patient_name)
                current_patient = self.Tree.AppendItem(root, s)

            if df.study_id != study_id:             # New study
                series_number = None
                study_id = df.study_id
                study_description = df.study_description
                s = "Study %s %s" % (df.study_id, df.study_description)
                current_study = self.Tree.AppendItem(current_patient, s)

            if df.series_number != series_number:   # New series
                series_number = df.series_number
                series_description = df.series_description
                s = "Series %d %s" % (df.series_number, df.series_description)
                current_series = self.Tree.AppendItem(current_study, s)

            description, icon_index = self.get_item_description(df)

            item_id = self.Tree.AppendItem(current_series,description,icon_index)
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

        
class DicomFileInfo(object):
    """
    A convenience container for a small subset of attributes from a DICOM file.

    The advantages of this over a pydicom Dataset object are --
    1) This uses Pythonic names
    2) All attributes are guaranteed to exist (although they may be empty)
    3) Allows flexible sorting. Lists of DicomFileInfo objects will sort on
       the attributes listed in the class-level attribute 'sort_attributes'.
       To sort on e.g. study_id and series_number:

           files = get_all_files(path)
           DicomFileInfo.sort_attributes = ('study_id', 'series_number')
           files.sort()

       One cannot sort on arbitrary dataset objects, only on instance attributes.

    """
    sort_attributes = ('patient_name','study_id','series_number','instance_number')

    def __init__(self, dataset, filename=""):
        """
        The dataset param must be a pydicom dataset. The filename is
        stored in this object for convenience but is otherwise unused.

        'slice_location' calculated from ImagePositionPatient and ImageOrientationPatient
        if present according to the algorithm outlined at great length by Jolinda Smith:
          http://www.itk.org/pipermail/insight-users/2003-September/004762.html
          http://www.cmake.org/pipermail/insight-users/2003-September/004762.html
        If those links die you can beseech Google:
          http://www.google.com/search?q=%22slightly-rotated+coronal+acquisition%22

        """
        self.dataset            = dataset
        self.filename           = filename
        self.patient_name       = str(dataset.PatientName)  # 0010 0010
        self.patient_id         = dataset.PatientID         # 0010 0020
        self.study_instance_uid = dataset.StudyInstanceUID  # 0020 000D
        self.study_id           = ''                        # 0020 0010
        self.study_description  = ''                        # 0008 1030
        self.image_id           = 0                         # 0054 0400
        self.instance_number    = 0                         # 0020 0013
        self.series_number      = 0                         # 0020 0011
        self.series_description = ''                        # 0008 103E
        self.series_instance_uid = ''                       # 0020 000E
        self.echo_time          = 0.001                     # 0018 0081
        self.acquisition_time   = 0                         # 0008 0032

        # Here I populate the not-always-present values
        if 'ImageType' in dataset:         self.image_type          = dataset.ImageType         # 0008 0008
        if 'StudyID' in dataset:           self.study_id            = dataset.StudyID
        if 'ImageID' in dataset:           self.image_id            = dataset.ImageID
        if 'SeriesNumber' in dataset:      self.series_number       = dataset.SeriesNumber
        if 'StudyDescription' in dataset:  self.study_description   = dataset.StudyDescription
        if 'SeriesInstanceUID' in dataset: self.series_instance_uid = dataset.SeriesInstanceUID
        if 'SeriesDescription' in dataset: self.series_description  = dataset.SeriesDescription
        if 'InstanceNumber' in dataset:    self.instance_number     = dataset.InstanceNumber
        if 'EchoTime' in dataset:          self.echo_time           = dataset.EchoTime
        if 'AcquisitionTime' in dataset:   self.acquisition_time    = dataset.AcquisitionTime

        self.slice_location = 0
        if  ('ImagePositionPatient'    in dataset)  and \
            ('ImageOrientationPatient' in dataset)  and \
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

        # I create a shortcut reference for the dataset
        ds = self.dataset

        if 'PixelData' not in ds:
            return None
        else:
            RAW_UIDS = (pydicom.uid.ExplicitVRLittleEndian,
                        pydicom.uid.ExplicitVRBigEndian)

            if ds.TransferSyntaxUID not in RAW_UIDS:
                raise NotImplementedError('Unsupported: TransferSyntaxUID (0002,0010) == %d' % ds.TransferSyntaxUID)
            else:
                # Data is in raw mode
                pixels = ds.pixel_array
                pixels = pixels.flatten().tolist()

            # We only handle monochrome images, so check PhotometricInterpretation
            # has Monochrome. If not in dataset, we assume MONOCHROME1 or MONOCHROME2.
            # MONOCHROME1 is inverse scale, so might have to convert. Finally
            # convert to RGB data string which wx.Image wants.

            if 'PhotometricInterpretation' in ds:
                if not ds.PhotometricInterpretation.startswith('MONOCHROME'):
                    msg = 'Unsupported: PhotometricInterpretation (0028,0004) == %d' % ds.PhotometricInterpretation
                    raise NotImplementedError(msg)
            else:
                if ('WindowCenter' in ds) and ('WindowWidth' in ds):
                    normalize_pixels(pixels, ds.WindowCenter, ds.WindowWidth)

                if 'PhotometricInterpretation' in ds:
                    if ds.PhotometricInterpretation == 'MONOCHROME1':
                        pixels = [_MAX_PIXEL - n for n in pixels]

                pixels = "".join([chr(int(round(n))) * 3 for n in pixels])
                image = wx.EmptyImage(ds.Rows, ds.Columns)
                image.SetData(pixels)

            return image



def normalize_pixels(pixels, window_center, window_width):
    """
    Applies a linear conversion to the pixels based on the Window
    Center (0028,1050) and Window Width (0028,1051) values of the dataset.
    Pixels must be an iterable (e.g. list); it is altered in-place.
    The pixels returned are floating point values between _MIN_PIXEL and
    _MAX_PIXEL.

    This implements the algorithm presented in the DICOM specification
    Part 3: Information Object Definitions, Annex C.11.2.1.2.
    ftp://medical.nema.org/medical/dicom/2008/08_03pu.pdf

    Here's the pseudocode from the DICOM spec:
        if      (x <= c - 0.5 - (w-1)/2), then y = ymin
        else if (x >  c - 0.5 + (w-1)/2), then y = ymax
        else y = ((x - (c - 0.5)) / (w-1) + 0.5) * (ymax - ymin) + ymin

    The Python code below is the same as the pseudocode with some loop
    invariants calculated here to improve performance.

    """
    min_threshold = (window_center - 0.5 - ((window_width - 1) / 2))
    max_threshold = (window_center - 0.5 + ((window_width - 1) / 2))
    output_range = _MAX_PIXEL - _MIN_PIXEL

    window_width  -= 1.0
    window_center -= 0.5

    for i, value in enumerate(pixels):
        if value <= min_threshold:
            pixels[i] = _MIN_PIXEL
        elif value > max_threshold:
            pixels[i] = _MAX_PIXEL
        else:
            pixels[i] = ((((value - window_center) / window_width) + 0.5) * output_range) + _MIN_PIXEL


def is_dicom(filename):
    """
    Returns True if the file in question is a DICOM file, else False. 
    
    - a DICOM file starts with 128 reserved bytes followed by 'DICM'.
    - ref: DICOM spec, Part 10: Media Storage and File Format for Media 
           Interchange, 7.1 DICOM FILE META INFORMATION 
    """
    if os.path.isfile(filename):
        f = open(filename, 'rb')
        s = f.read(132)
        if isinstance(s,(bytes, bytearray)):
            try:
                s = s.decode('utf-8')
            except:
                # bugfix - was trying to read PNG and found a char utf8 did not like
                try:
                    s = s.decode('utf-16')
                except:
                    f.close()
                    return False
        f.close()
        return s.endswith('DICM')
    else:
        return False


def get_all_files(path):
    """
    Gets all DICOM files in the directory indicated by the path and
    returns them as an unsorted list of DicomFileInfo instances.
    To sort the files, call the sort() method of the returned list.
    """
    return [f for f in get_files(path)]


def get_files(path):
    """
    Return DicomFileInfo instances (or a subclass thereof) for the DICOM files
    in the directory indicated by the path. This is a Python generator so they
    are returned one by one.
    
    Example, use it in a loop like so:

        for dicom_file in get_files(the_path):
            do_something(dicom_file)

    See also get_all_files().

    """
    # list files in path, add path to fname, remove if not a file
    filenames = os.listdir(path)
    filenames = [os.path.join(path, filename) for filename in filenames]
    filenames = [filename for filename in filenames if os.path.isfile(filename)]

    for filename in filenames:
        if is_dicom(filename):
            dataset = pydicom.dcmread(filename)
            dataset.decode()                        # change strings to unicode 
            df = DicomFileInfo(dataset, filename)
            yield df
        #else:
            # Not a DICOM file; ignore it.        


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


# =============================================================================    
# Dialog DICOM Browser Subclasses that show Manufacturer Specific Files

class SiemensSpectroscopyDicomBrowserDialog(DicomBrowserDialog):
    def __init__(self, parent=None, multi_select=False, path=''):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        self.manufacturer = ['SIEMENS', ]
        DicomBrowserDialog.__init__(self, parent, None, multi_select, path=path)


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
        
        if TAG_SOP_CLASS_UID in dataset:
            item_uid = dataset[TAG_SOP_CLASS_UID].value.upper() 
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
    


class PhilipsSpectroscopyDicomBrowserDialog(DicomBrowserDialog):
    def __init__(self, parent=None, multi_select=False, path=''):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        self.manufacturer=['PHILIPS',]
        DicomBrowserDialog.__init__(self, parent, None, multi_select, path=path)


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
        
        return util_philips.is_mrs_dicom(dicom_file.dataset)      
#!/usr/bin/env python
"""
Expansion of matplotlib embed in wx example  
see http://matplotlib.org/examples/user_interfaces/embedding_in_wx4.html

This version, image_panel_toolbar.py, is an expansion that allows users to 
display one or more images in vertical arrangement of axes. A toolbar with
custom icons is attached to the bottom that allow the user to pan, zoom,
save figures, and undo/home the images in the axes. Note that standard
settings are for multiple axes to be linked for pan/zoom operations.  

Brian J. Soher, Duke University, April, 2014

Update - 31 July 2022 - Py3 and Newer wx

- all .lines = [] to .lines.clear()
- all .patches = [] to .patches.clear()
- all .images = [] to .images.clear()


"""
# Python modules

import os

# third party modules
import wx
import numpy as np
import matplotlib
import matplotlib.cm as cm

matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backend_bases          import NavigationToolbar2, cursors
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from matplotlib.figure    import Figure
from numpy.random         import rand
from wx.lib.embeddedimage import PyEmbeddedImage


LEVMAX =  383
LEVMIN = -127
WIDMAX =  255
WIDMIN =    0
LEVSTR =  128


# MPL widlev example
#
#http://matplotlib.1069221.n5.nabble.com/Pan-zoom-interferes-with-user-specified-right-button-callback-td12186.html


class ImagePanelToolbar2(wx.Panel):
    """
    The ImagePanel has a Figure and a Canvas and 'n' Axes. The user defines
    the number of axes on Init and this number cannot be changed thereafter.
    However, the user can change the number of axes displayed in the Figure.
    
    Axes are specified on Init because the zoom and widlev actions 
    need an axes to attach to initialize properly.  
    
    on_size events simply set a flag, and the actual resizing of the figure is 
    triggered by an Idle event.
    
    """

    # Set _EVENT_DEBUG to True to activate printing of messages to stdout 
    # during events.
    _EVENT_DEBUG = False

    def __init__(self, parent, naxes=1, 
                               data=None, 
                               cmap=cm.gray, 
                               color=None,
                               bgcolor="#ffffff",
                               vertOn=False,
                               horizOn=False,
                               lcolor='gold', 
                               lw=0.5,
                               layout='vertical',
                               **kwargs):
        
        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        self.figure = Figure(figsize=(4,4), dpi=100)

        self.cmap       = [cm.gray for i in range(naxes)]
        self.imageid    = [None    for i in range(naxes)]
        self.width      = [WIDMAX  for i in range(naxes)]
        self.level      = [LEVSTR  for i in range(naxes)]
        self.vmax       = [WIDMAX  for i in range(naxes)]
        self.vmin       = [LEVSTR  for i in range(naxes)]
        self.vmax_orig  = [WIDMAX  for i in range(naxes)]
        self.vmin_orig  = [LEVSTR  for i in range(naxes)]
        
        # here we create the required naxes, add them to the figure, but we
        # also keep a permanent reference to each axes so they can be added
        # or removed from the figure as the user requests 1-N axes be displayed
        self.naxes = naxes   
        self.axes  = []

        if layout == 'vertical':
        self.axes.append(self.figure.add_subplot(naxes,1,1))
        if naxes > 1:
                for i in range(naxes-1):
                    self.axes.append(self.figure.add_subplot(naxes,1,i+2, sharex=self.axes[0], sharey=self.axes[0]))
        else:
            self.axes.append(self.figure.add_subplot(1,naxes,1))
            if naxes > 1:
                for i in range(naxes - 1):
                    self.axes.append(self.figure.add_subplot(1,naxes,i+2, sharex=self.axes[0], sharey=self.axes[0]))

        self.all_axes = list(self.axes)
 
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.figure.subplots_adjust(left=0.0,right=1.0,
                                    bottom=0.0,top=1.0,
                                    wspace=0.0,hspace=0.0)

        # for images we don't show x or y axis values
        for axes in self.axes:
            axes.set_adjustable('box')
            axes.xaxis.set_visible(False)
            axes.yaxis.set_visible(False)

        self.set_color( color )
        self.set_axes_color()

        if not data or len(data) != naxes:
            data = self._default_data()
        self.set_data(data)
        self.update(no_draw=True)   # we don't have a canvas yet

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        if layout == 'vertical':
            self.sizer_canvas = wx.BoxSizer(wx.VERTICAL)
            self.sizer_canvas.Add(self.canvas, 1, wx.EXPAND, 10)
        else:
            self.sizer_canvas = wx.BoxSizer(wx.HORIZONTAL)
            self.sizer_canvas.Add(self.canvas, 1, wx.EXPAND, 10)

        self.sizer.Add(self.sizer_canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)

        # Capture the paint message
#        wx.EVT_PAINT(self, self.on_paint)
        wx.EvtHandler.Bind(self, wx.EVT_PAINT, self.on_paint)

        self.toolbar = NavigationToolbar3Wx(self.canvas, 
                                            self, 
                                            vertOn=vertOn, 
                                            horizOn=horizOn,
                                            lcolor=lcolor,
                                            lw=lw)

        #self.toolbar = NavigationToolbar2Wx(self.canvas)

        self.toolbar.Realize()
        if wx.Platform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
#            tw, th = self.toolbar.GetSize()
#            fw, fh = self.canvas.GetSize()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
#            self.toolbar.SetSize(wx.Size(fw, th))
            #self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
            #self.sizer.Add(self.toolbar, 0, wx.ALIGN_CENTER | wx.EXPAND)
            self.sizer.Add(self.toolbar, 0, wx.EXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()
        
        self.SetSizer(self.sizer)
        self.Fit()

        # Always link scroll events, but only handle them if on_scroll()
        # is overloaded by the user.
        self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)

        self.keypress_id   = self.canvas.mpl_connect('key_press_event',   self.on_key_press)
        self.keyrelease_id = self.canvas.mpl_connect('key_release_event', self.on_key_release)
        self.shift_is_held = False


    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================

    def on_key_press(self, event):
        if event.key == 'shift':
            self.shift_is_held = True
    
    def on_key_release(self, event):
        if event.key == 'shift':
            self.shift_is_held = False

    def _dprint(self, a_string):
        if self._EVENT_DEBUG:
            print(a_string)

    def _on_scroll(self, event):
        """
        This is the internal method that organizes the data that is sent to the
        external user defined event handler for scroll events. In here we 
        determine which axis we are in, then call the (hopefully) overloaded 
        on_scroll() method
        
        """
        if event.inaxes == None: return
        for i,axes in enumerate(self.axes):
            if axes == event.inaxes:
                iplot = i
        self.on_scroll(event.button, event.step, iplot)        


    def _default_data(self):
        data = []
        for i in range(self.naxes):
            data.append([self._dist(128),])
        return data


    def _dist(self, n, m=None):  
        """
        Return a rectangular array in which each pixel = euclidian
        distance from the origin.

        """
        n1 = n

        if not m:
            m1 = n
        else:
            m1 = m

        x = np.arange(n1)
        x = np.array([val**2 if val < (n1-val) else (n1-val)**2 for val in x ])

        a = np.ndarray((n1,m1),float)   #Make array

        for i in range(int((m1/2)+1)):       #Row loop
            y = np.sqrt(x + i**2.0)      #Euclidian distance
            a[i,:] = y              #Insert the row
            if i != 0:
                a[m1-i,:] = y      #Symmetrical

        return a
        

    def on_paint(self, event):
        # this is necessary or the embedded MPL canvas does not show
        self.canvas.draw()
        event.Skip()


    #=======================================================
    #
    #           Default Event Handlers  
    #
    #=======================================================
        
    def on_scroll(self, button, step, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_scroll,     button='+str(button)+'  step='+str(step)+'  Index = '+str(iplot))

    def on_motion(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_motion,          xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))
        
    def on_select(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_select,          xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

    def on_panzoom_release(self, xloc, yloc):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_panzoom_release, xloc='+str(xloc)+'  yloc='+str(yloc))
        
    def on_panzoom_motion(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_panzoom_motion,  xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

    def on_level_press(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_level_press,     xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

    def on_level_release(self, xloc, yloc):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_level_release,   xloc='+str(xloc)+'  yloc='+str(yloc))

    def on_level_motion(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_level_motion,    xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

    #=======================================================
    #
    #           User Accessible Data Functions  
    #
    #=======================================================

    def set_data(self, data, index=None, keep_norm=False):
        """
        User can set data into one or all axes using this method.
        
        If index is supplied, we assume that only one ndarray is being
        passed in via the data parameter. If no index is supplied then we
        assume that a list of ndarrays the size of self.naxes is being 
        passed in to replace all data in all axes.

        Example 1 - Data is a list of dicts
        
            raw  = {'data'  : raw_data,        # 2D numpy array 
                    'alpha' : 1.0 }            # value 0.0-1.0

            fit  = {'data' : fit_data,         # 2D numpy array 
                    'alpha' : 1.0 }            # value 0.0-1.0

            data = [raw, fit]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized, no_draw=True)
            self.view.canvas.draw()

        Example 2 - Data is a single numpy array, the colors dict will use
                    default values set in set_data() method
        
            data = [raw_data,]          # 2D numpy array
            data = [[data]]
            self.view.set_data(data)    # alpha defaults to 1.0
            self.view.update(set_scale=not self._scale_intialized, no_draw=True)
            self.view.canvas.draw()
                    
        """
        for i, item in enumerate(data):
            for j, dat in enumerate(item):
                if isinstance(dat, dict):
                    # Dict in this item, but ensure all keys are present
                    if 'data' not in list(dat.keys()):
                        raise ValueError("must have a data array in the dict sent to set_data()")
                    if 'alpha' not in list(dat.keys()):
                        dat['alpha'] = 1.0
                    if 'cmap' not in list(dat.keys()):
                        dat['cmap'] = self.cmap[i]
                    if 'vmax' not in list(dat.keys()):
                        dat['vmax'] = dat['data'].max()
                    if 'vmin' not in list(dat.keys()):
                        dat['vmin'] = dat['data'].min()
                    if 'keep_norm' not in list(dat.keys()):
                        dat['keep_norm'] = keep_norm
                    if 'patches' not in list(dat.keys()):
                        dat['patches'] = None
                    if 'lines' not in dat.keys():
                        dat['lines'] = None
                else:
                    # Only data in this item, so add all default values 
                    dat = { 'data'      : dat,
                            'alpha'     : 1.0,
                            'cmap'      : self.cmap[i],
                            'vmax'      : dat.max(),
                            'vmin'      : dat.min(),
                            'keep_norm' : keep_norm,
                            'patches'   : None,
                            'lines'     : None
                          }
                        
                item[j] = dat
        
        if index:
            if index < 0 or index >= self.naxes:
                raise ValueError("index must be within that number of axes in the plot_panel")
            
            if data[0][0]['data'].shape != self.data[0][0]['data'].shape:
                raise ValueError("new data must be a same number of spectral points as existing data")
            
            # even though we are inserting into an index, I want to force users
            # to submit a dict in a list of lists format so it is consistent 
            # with submitting a whole new set of data (below). We just take the
            # first list of dicts from the submitted data and put it in the 
            # index position
            self.data[index] = data[0]
            
        else:
            if len(data) != self.naxes:
                raise ValueError("data must be a list with naxes number of ndarrays")
            for item in data:
                for dat in item:
                    d = dat['data']
                
                    padding = 2 - len(d.shape)
                    if padding > 0:
                        d.shape = ([1] * padding) + list(d.shape)
                    elif padding == 0:
                        # Nothing to do; data already has correct number of dims
                        pass
                    else:
                        # padding < 0 ==> data has too many dims
                        raise ValueError("Data with shape %s has too many dimensions" % str(item.shape))
                
                    if d.shape != data[0][0]['data'].shape:
                        raise ValueError("all ndarrays must have same dimensions")
            
            self.data = data

        for i,data in enumerate(self.data):
            if not data[0]['keep_norm']:
                self.vmax_orig[i] = data[0]['vmax']
                self.vmin_orig[i] = data[0]['vmin']
                self.vmax[i]      = data[0]['vmax']
                self.vmin[i]      = data[0]['vmin']
                
#         if not keep_norm:
#             for i,data in enumerate(self.data):
#                 self.vmax_orig[i] = data[0]['vmax']
#                 self.vmin_orig[i] = data[0]['vmin']
#                 self.vmax[i]      = data[0]['vmax']
#                 self.vmin[i]      = data[0]['vmin']


    def update(self, index=None, keep_norm=False, no_draw=False):
        """
        Convenience function that runs through all the typical steps needed
        to refresh the screen after a set_data().
        
        The set_scale option is typically used only once to start set the 
        bounding box to reasonable bounds for when a zoom box zooms back 
        out.
        
        """
        self.update_images(index=index)
        if not no_draw:
            self.canvas.draw()

    
    def update_images(self, index=None, force_bounds=False):
        """
        Sets the data from the normalized image numpy arrays into the axes.
        
        We also set the axes dataLim and viewLim ranges here just in case
        the dimensions on the image being displayed has changed. This reset
        allows the zoom reset to work properly.
        
        """
        indices = self.parse_indices(index)
        
        for i in indices:

            axes = self.all_axes[i]
            
            if axes.images:
                yold, xold = axes.images[0].get_array().shape
            else:
                yold, xold = -1,-1
            
            axes.images.clear()
            
            ddict    = self.data[i][0]
            data     = ddict['data'].copy() 
            alpha    = ddict['alpha']
            cmap     = ddict['cmap']
            vmax     = self.vmax[i]
            vmin     = self.vmin[i]
            patches  = ddict['patches']
            lines    = ddict['lines']
                
            xmin, xwid, ymin, ywid = 0, data.shape[1], 0, data.shape[0]
            
            # Explicit set origin here so the user rcParams value is not used. 
            # This keeps us consistent across users.
            self.imageid[i] = axes.imshow(data, cmap=cmap, 
                                                alpha=alpha, 
                                                vmax=vmax, 
                                                vmin=vmin, 
                                                aspect='equal', 
                                                origin='upper') 

            if patches is not None:
                self.imageid[i].axes.patches = []
                for patch in patches:
                    self.imageid[i].axes.add_patch(patch)
            else:
                self.imageid[i].axes.patches.clear()
            
            if len(self.imageid[i].axes.lines) > 2:
                # should be two lines in here for cursor tracking
                self.imageid[i].axes.lines = self.imageid[i].axes.lines[0:2]
            if lines is not None:
                for line in lines:
                    self.imageid[i].axes.add_line(line)

            if xold != xwid or yold!=ywid or force_bounds: 
                xmin -= 0.5     # this centers the image range to voxel centers
                ymin -= 0.5   
                # Set new bounds for dataLims to x,y extent of data in the 
                # new image. On reset zoom this is how far we reset the limits.
                axes.ignore_existing_data_limits = True
                axes.update_datalim([[xmin,ymin],[xmin+xwid,ymin+ywid]])
#                 # Matches viewLims view limits to the new data. By 
#                 # default, new data and updating causes display to show the 
#                 # entire image. Any zoom setting is lost.
#                 axes.set_xlim((xmin, xmin+xwid), auto=None)
#                 axes.set_ylim((ymin, ymin+ywid), auto=None)
# # 
# #             # may need this line in future if we do overlay
# #                 self.figure.hold(True)
# #             self.figure.hold(False)


    def parse_indices(self, index=None):
        """ 
        Ensure we know what data axes to act upon
         - index can be a list or scalar
         - if list, must be naxes items or less
         - list/scalar values need to be in range of 0 to naxes-1
        
        """
        if index is None:
            indices = list(range(self.naxes))
        else:
            if isinstance(index, list):
                if len(index) <= self.naxes:
                    if all(index < self.naxes):
                        indices = index
                    else:
                        raise ValueError("index in list outside naxes range")
                else:
                    raise ValueError("too many index entries")
            else:
                if index < self.naxes:
                    indices = [index]
                else:
                    raise ValueError("scalar index outside naxes range")
            
        return indices 


    #=======================================================
    #
    #           User Accessible Display Functions  
    #
    #=======================================================

    def set_color( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )
        
    def set_axes_color( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        for axis in self.all_axes:
            axis.spines['bottom'].set_color(clr)
            axis.spines['top'].set_color(clr) 
            axis.spines['right'].set_color(clr)
            axis.spines['left'].set_color(clr)        








#------------------------------------------------------------------------------
# NavigationToolbar3Wx
#
# I did not inherit from NavigationToolbar2Wx because that version is only
# minimally different from NavigationToolbar2. And I wanted to add leveling
# events which required access to many of the methods.  Thus, monkey-patching
# is the way to go.
#
# This toolbar is specific to use in canvases that display images.
#
#------------------------------------------------------------------------------


class NavigationToolbar3Wx(NavigationToolbar2, wx.ToolBar):

    # list of toolitems to add to the toolbar, format is:
    # (
    #   text, # the text of the button (often not visible to users)
    #   tooltip_text, # the tooltip shown on hover (where possible)
    #   image_file, # name of the image for the button (without the extension)
    #   name_of_method, # name of the method in NavigationToolbar2 to call
    # )
    toolitems = (
        ('Home', 'Reset original view', 'nav3_home', 'home'),
        ('Back', 'Back to  previous view', 'nav3_back', 'back'),
        ('Forward', 'Forward to next view', 'nav3_forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'nav3_move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'nav3_zoom_to_rect', 'zoom'),
        ('Level','Set width/level with right mouse', 'nav3_contrast', 'level'),
        ('Cursors', 'Track mouse movement with crossed cursors', 'nav3_crosshair', 'cursors'),
        (None, None, None, None),
        ('Subplots', 'Configure subplots', 'nav3_subplots', 'configure_subplots'),
        ('Save', 'Save the figure', 'nav3_filesave', 'save_figure'),
      )
    
    
    def __init__(self, canvas, parent, vertOn=False, horizOn=False, lcolor='gold', lw=0.5, coordinates=True):
        wx.ToolBar.__init__(self, canvas.GetParent(), -1)

        if 'wxMac' in wx.PlatformInfo:
            self.SetToolBitmapSize((24, 24))
        self.wx_ids = {}
        for text, tooltip_text, image_file, callback in self.toolitems:

            if text is None:
                self.AddSeparator()
                continue

            bmp = nav3_catalog[image_file].GetBitmap()

            self.wx_ids[text] = wx.NewIdRef()
            if text in ['Pan', 'Level', 'Cursors']:
                self.AddCheckTool(self.wx_ids[text], ' ', bmp, shortHelp=text, longHelp=tooltip_text)
            elif text in ['Zoom', 'Subplots']:
                pass  # don't want this in my toolbar
            else:
                self.AddTool(self.wx_ids[text], ' ', bmp, text)

            self.Bind(wx.EVT_TOOL, getattr(self, callback), id=self.wx_ids[text])

        self._coordinates = coordinates
        if self._coordinates:
            self.AddStretchableSpace()
            self._label_text = wx.StaticText(self)
            self.AddControl(self._label_text)

        self.Realize()

        NavigationToolbar2.__init__(self, canvas)

        # bjs-start

        self.canvas = canvas
        self.parent = parent
        self._idle = True
        self.statbar = self.parent.statusbar

        self._button_pressed = None
        self._xypress = None

        # turn off crosshair cursors when mouse outside canvas
        self._idAxLeave  = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self._idFigLeave = self.canvas.mpl_connect('figure_leave_event', self.leave)
        self._idRelease = self.canvas.mpl_connect('button_release_event', self.release_local)
        self._idPress = None
        self._idDrag = self.canvas.mpl_connect( 'motion_notify_event', self.mouse_move)
        
        # set up control params for cursor crosshair functionality
        self._cursors = False
        self.vertOn  = vertOn
        self.horizOn = horizOn
        self.vlines = []
        self.hlines = []
        if vertOn:  self.vlines = [ax.axvline(0, visible=False, color=lcolor, lw=lw) for ax in self.parent.axes]
        if horizOn: self.hlines = [ax.axhline(0, visible=False, color=lcolor, lw=lw) for ax in self.parent.axes]

        
    def get_canvas(self, frame, fig):
        # saw this in NavigationToolbar2WxAgg, so included here
        return FigureCanvasWxAgg(frame, -1, fig)


    def zoom(self, *args):
            
        for item in ['Pan', 'Level']:
            if item in list(self.wx_ids.keys()): self.ToggleTool(self.wx_ids[item], False)
            
        # if 'Pan' in list(self.wx_ids.keys()):
        #     self.ToggleTool(self.wx_ids['Pan'], False)
        # if 'Level' in list(self.wx_ids.keys()):
        #     self.ToggleTool(self.wx_ids['Level'], False)

        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
        NavigationToolbar2.zoom(self, *args)
        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_connect('button_release_event', self.release_local)

 
    def pan(self, *args):

        for item in ['Zoom', 'Level']:
            if item in list(self.wx_ids.keys()): self.ToggleTool(self.wx_ids[item], False)

        # if 'Zoom' in list(self.wx_ids.keys()):
        #     self.ToggleTool(self.wx_ids['Zoom'], False)
        # if 'Level' in list(self.wx_ids.keys()):
        #     self.ToggleTool(self.wx_ids['Level'], False)

        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
        NavigationToolbar2.pan(self, *args)
        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_connect('button_release_event', self.release_local)


    def drag_pan(self, event):
        NavigationToolbar2.drag_pan(self, event)
        self.mouse_move(event)
        
        
    def level(self, *args):
        """Activate the width/level tool. change values with right button"""
        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)

        
        if self.parent.axes[0].get_navigate_mode() == 'LEVEL':
            if 'Level' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Level'], False)
        else: 
            if 'Level' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Level'], True)
            if 'Pan' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Pan'], False)
            if 'Zoom' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Zoom'], False)      
        
        # set the pointer icon and button press funcs to the
        # appropriate callbacks
    
        if self.parent.axes[0].get_navigate_mode() == 'LEVEL':
            self.parent.axes[0].set_navigate_mode(None)
        else:
            self.parent.axes[0].set_navigate_mode('LEVEL')

        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''
    
        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''
    
        if self.parent.axes[0].get_navigate_mode() == 'LEVEL':
            self._idPress   = self.canvas.mpl_connect( 'button_press_event',   self.press_level)
            self._idRelease = self.canvas.mpl_connect( 'button_release_event', self.release_level)
            self.mode = 'width/level'
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)
    
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self.parent.axes[0].get_navigate_mode())
    
        self.set_message(self.mode)
 
        if self.parent.axes[0].get_navigate_mode() is None:
            self._idRelease = self.canvas.mpl_connect('button_release_event', self.release_local)

 
    def press_level(self, event):
        """the press mouse button in width/level mode callback"""
        
        if event.button == 3:
            self._button_pressed = 3
        else:
            self._button_pressed = None
            return
        
        x, y = event.x, event.y
        
        # push the current view to define home if stack is empty
        if self._nav_stack.empty():
            self.push_current()
        
        self._xypress = []
        for i, a in enumerate(self.canvas.figure.get_axes()):
            if (x is not None and y is not None and a.in_axes(event) and a.get_navigate() and a.can_pan()):
                self._xypress.append([a, i, event.x, event.y])
                self.canvas.mpl_disconnect(self._idDrag)
                self._idDrag = self.canvas.mpl_connect('motion_notify_event', self.drag_level)
        
        self.press_local(event)
 
 
    def drag_level(self, event):
        """the drag callback in width/level mode"""
        
        for a, indx, xprev, yprev in self._xypress:
            
            xdelt = int((event.x-xprev))   
            ydelt = int((event.y-yprev))  
        
            if abs(ydelt) >= abs(xdelt):    
                self.parent.level[indx] += ydelt
            else:                           
                self.parent.width[indx] += xdelt
    
            vmax0 = self.parent.vmax_orig[indx] 
            vmin0 = self.parent.vmin_orig[indx]
            wid   = max(WIDMIN, min(WIDMAX, self.parent.width[indx]))
            lev   = max(LEVMIN, min(LEVMAX, self.parent.level[indx]))
            vtot  = vmax0 - vmin0
            vmid  = vmin0 + vtot * (lev/WIDMAX)
            vwid  = vtot * (wid/WIDMAX)
            vmin = vmid - (vwid/2.0)
            vmax = vmid + (vwid/2.0)
            
            self.parent.vmax[indx]  = vmax
            self.parent.vmin[indx]  = vmin
            self.parent.width[indx] = wid   # need this in case values were
            self.parent.level[indx] = lev   # clipped by MIN MAX bounds
            
            self.parent.update_images(index=indx)

            # prep new 'last' location for next event call
            self._xypress[0][2] = event.x
            self._xypress[0][3] = event.y
            
        self.mouse_move(event)
    
        self.dynamic_update()

 
    def release_level(self, event):
        """the release mouse button callback in width/level mode"""
        
        if self._button_pressed is None:
            return
        self.canvas.mpl_disconnect(self._idDrag)
        self._idDrag = self.canvas.mpl_connect( 'motion_notify_event', self.mouse_move)
        for a, ind, xlast, ylast in self._xypress:
            pass
        if not self._xypress:
            return
        self._xypress = []
        self._button_pressed = None
        self.push_current()
        self.release_local(event)
        self.canvas.draw()


    def cursors(self, *args):
        """
        Toggle the crosshair cursors tool to show vertical and 
        horizontal lines that track the mouse motion, or not.
        
        """
        if self._cursors:
            self._cursors = False
            if 'Cursors' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Cursors'], False)
        else: 
            self._cursors = True
            if 'Cursors' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Cursors'], True)

 
    def configure_subplots(self, evt):
        frame = wx.Frame(None, -1, "Configure subplots")
    
        toolfig = Figure((6,3))
        canvas = self.get_canvas(frame, toolfig)
    
        # Create a figure manager to manage things
        figmgr = FigureManager(canvas, 1, frame)
    
        # Now put all into a sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
        # This way of adding to sizer allows resizing
        sizer.Add(canvas, 1, wx.LEFT|wx.TOP|wx.GROW)
        frame.SetSizer(sizer)
        frame.Fit()
        tool = SubplotTool(self.canvas.figure, toolfig)
        frame.Show()
 
    def save_figure(self, *args):
        # Fetch the required filename and file type.
        filetypes, exts, filter_index = self.canvas._get_imagesave_wildcards()
        default_file = self.canvas.get_default_filename()
        dlg = wx.FileDialog(self._parent, "Save to file", "", default_file,
                            filetypes,
                            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        dlg.SetFilterIndex(filter_index)
        if dlg.ShowModal() == wx.ID_OK:
            dirname  = dlg.GetDirectory()
            filename = dlg.GetFilename()
            format = exts[dlg.GetFilterIndex()]
            basename, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if ext in ('svg', 'pdf', 'ps', 'eps', 'png') and format!=ext:
                #looks like they forgot to set the image type drop
                #down, going with the extension.
                warnings.warn('extension %s did not match the selected image type %s; going with %s'%(ext, format, ext), stacklevel=0)
                format = ext
            try:
                self.canvas.print_figure(
                    os.path.join(dirname, filename), format=format)
            except Exception as e:
                error_msg_wx(str(e))
 
    def set_cursor(self, cursor):
        cursor =wx.Cursor(cursord[cursor])
        self.canvas.SetCursor( cursor )

    def press_local(self, event):
        
        xloc, yloc = self.get_bounded_xyloc(event)        
        if len(self._xypress)>0:
        item = self._xypress[0]
        axes, iplot = item[0], item[1]
        if self.mode == 'pan/zoom':
            for line in self.vlines + self.hlines:
                line.set_visible(False)
            self.dynamic_update()
            pass
        elif self.mode == 'zoom rect':
            for line in self.vlines + self.hlines:
                line.set_visible(False)
            self.dynamic_update()
            pass
        elif self.mode == 'width/level':
            for line in self.vlines + self.hlines:
                line.set_visible(False)
            self.dynamic_update()
            self.parent.on_level_press(xloc, yloc, iplot)
        elif self.mode == '':
            # no toggle buttons on
            # but, maybe we want to show crosshair cursors
            if self._cursors:
                for line in self.vlines + self.hlines:
                    line.set_visible(True)
            self.dynamic_update()
            pass
        else:
            # catch all
            pass
    
    def release_local(self, event):
        # legacy code from NavigationToolbar2Wx
        try: 
            del self.lastrect
        except AttributeError: 
            pass

        # find out what plot we released in
        iplot = None
        for i,axes in enumerate(self.parent.axes):
            if axes == event.inaxes:
                iplot = i
        
        # get bounded location and call user event
        xloc, yloc = self.get_bounded_xyloc(event)
        if self.mode == 'pan/zoom':
            self.parent.on_panzoom_release(xloc, yloc)
        elif self.mode == 'zoom rect':
            pass
        elif self.mode == 'width/level':
            self.parent.on_level_release(xloc, yloc)
        elif self.mode == '':
            # no toggle buttons on
            self.parent.on_select(xloc, yloc, iplot)
        else:
            # catch all
            self.parent.on_select(xloc, yloc)


    def leave(self, event):
        """ Turn off the cursors as we move mouse outside axes or figure """

        for line in self.vlines+self.hlines: line.set_visible(False)
        self.dynamic_update()            
            

    def on_motion(self, event):
        
        # The following limit when motion events are triggered 
        
        if not event.inaxes and self.mode == '': 
            return
        if not event.inaxes and not self._button_pressed:
            # When we are in pan or zoom or level and the mouse strays outside
            # the canvas, we still want to have events even though the xyloc
            # can not be properly calculated. We are reporting other things 
            # that do not need xylocs for like width/level values. 
            return
        
        xloc, yloc = self.get_bounded_xyloc(event)

        iplot = None              
        if self._xypress:
            item = self._xypress[0]
            axis, iplot = item[0], item[1] 
        else:
            axis = None
            for i,axes in enumerate(self.parent.axes):
                if axes == event.inaxes:
                    iplot = i
        
        if self.mode == 'pan/zoom' and (self._button_pressed == 1 or self._button_pressed == 3):
            if iplot is not None:
                self.parent.on_panzoom_motion(xloc, yloc, iplot)
        elif self.mode == 'zoom rect' and (self._button_pressed == 1 or self._button_pressed == 3):
            if iplot is None:
                pass
        elif self.mode == 'width/level' and self._button_pressed == 3:
            if iplot is not None:
                self.parent.on_level_motion(xloc, yloc, iplot)
        elif self.mode == '':
            if iplot is not None:
                # no toggle buttons on
                self.parent.on_motion(xloc, yloc, iplot)

                if self._cursors and self.vertOn:
                    for line in self.vlines:
                        line.set_xdata((xloc, xloc))
                        line.set_visible(True)
                if self._cursors and self.horizOn:
                    for line in self.hlines:
                        line.set_ydata((yloc, yloc))
                        line.set_visible(True)    
                self.dynamic_update()            
        else:
            if iplot is not None:
                # catch all
                self.parent.on_motion(xloc, yloc, iplot)
            
  
    def mouse_move(self, event):
        # call base event
        NavigationToolbar2.mouse_move(self, event)

        if event.inaxes and self.parent.axes[0].get_navigate_mode():
            if (self.parent.axes[0].get_navigate_mode() == 'LEVEL' and self._lastCursor != 4):
                self.set_cursor(4)
                self._lastCursor = 4
        
        # call additional event
        self.on_motion(event)


    def get_bounded_xyloc(self, event):
        if not event.inaxes:
            return 0,0

        xloc, yloc = event.xdata, event.ydata
        
        # bound these values to be inside the size of the image.
        # - a pan event could yield negative locations.
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
        xmin,xmax = x0+0.5, x0+0.5+x1-1
        ymin,ymax = y0+0.5, y0+0.5+y1-1     # position swap due to 0,0 location 'upper'
        xloc = max(0, min(xmax, xloc))
        yloc = max(0, min(ymax, yloc))

        return xloc,yloc
    
    
    def dynamic_update(self):
        d = self._idle
        self._idle = False
        if d:
            self.canvas.draw()
            self._idle = True

    def draw_rubberband(self, event, x0, y0, x1, y1):
        'adapted from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/189744'
        canvas = self.canvas
        dc = wx.ClientDC(canvas)
    
        # Set logical function to XOR for rubberbanding
        dc.SetLogicalFunction(wx.XOR)
    
        # Set dc brush and pen
        # Here I set brush and pen to white and grey respectively
        # You can set it to your own choices
    
        # The brush setting is not really needed since we
        # dont do any filling of the dc. It is set just for
        # the sake of completion.
    
        wbrush =wx.Brush(wx.Colour(255,255,255), wx.TRANSPARENT)
        wpen =wx.Pen(wx.Colour(200, 200, 200), 1, wx.SOLID)
        dc.SetBrush(wbrush)
        dc.SetPen(wpen)
    
        dc.ResetBoundingBox()
        dc.BeginDrawing()
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
    
        if y1<y0: y0, y1 = y1, y0
        if x1<y0: x0, x1 = x1, x0
    
        w = x1 - x0
        h = y1 - y0
    
        rect = int(x0), int(y0), int(w), int(h)
        try: lastrect = self.lastrect
        except AttributeError: pass
        else: dc.DrawRectangle(*lastrect)  #erase last
        self.lastrect = rect
        dc.DrawRectangle(*rect)
        dc.EndDrawing()
    
    def set_status_bar(self, statbar):
        self.statbar = statbar
    
    def set_message(self, s):
       if self.statbar is not None:
           pass 
           #self.statbar.SetStatusText(s,0)
    
    def set_history_buttons(self):
        can_backward = (self._nav_stack._pos > 0)
        can_forward = (self._nav_stack._pos < len(self._nav_stack._elements) - 1)
        self.EnableTool(self.wx_ids['Back'], can_backward)
        self.EnableTool(self.wx_ids['Forward'], can_forward)
        # try:
        #     self.EnableTool(self.wx_ids['Back'], can_backward)
        #     self.EnableTool(self.wx_ids['Forward'], can_forward)
        # except:
        #     bob = 10

# ***************** Shortcuts to wxPython standard cursors************

cursord = { cursors.MOVE          : wx.CURSOR_HAND,
            cursors.HAND          : wx.CURSOR_HAND,
            cursors.POINTER       : wx.CURSOR_ARROW,
            cursors.SELECT_REGION : wx.CURSOR_CROSS,
            4                     : wx.CURSOR_BULLSEYE,
          }        


# ***************** Embed Icon Catalog Starts Here *******************

nav3_catalog = {}
nav3_index = []


#----------------------------------------------------------------------
nav3_crosshair = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAIBJREFUWEftlVEKgCAQBfcgXau/TtB9OkF/XbJ8sBsimlToFryB"
    "+RJxQGXl70yqG5vqBgMYwAAGrGpXhuBYEGvNwUF7Qaw1xwJSXgVgotmDuhL3XQtYgum+nHPw"
    "xD3gDrWA5lhAzi4B7t8wBvcN3bAH5QYDGMAABmCiPZ5qH0DkAEnpVB1zVLSlAAAAAElFTkSu"
    "QmCC")
nav3_index.append('nav3_crosshair')
nav3_catalog['nav3_crosshair'] = nav3_crosshair
getnav3_crosshairData = nav3_crosshair.GetData
getnav3_crosshairImage = nav3_crosshair.GetImage
getnav3_crosshairBitmap = nav3_crosshair.GetBitmap
getnav3_crosshairIcon = nav3_crosshair.GetIcon


#----------------------------------------------------------------------
nav3_back = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAmNJREFUSEutlbuL1FAUxhPB3dUBFRlmnC1URIvNYx55zmtzkgXZ"
    "cWYndq4gCGJtY+VfYCNoZ2FhY2lloY0gdlY2i1hspwhi4WPEKVZ30e9mb+KdMUrGzcCPJOee"
    "831z7wknkuGFBdOsy4ZRk9nVcYzoGmNZDdnwhgWTQhP0warphaqxPNgf1zUataTOsn7XSuxn"
    "0vl1i4bXfL81z1lg1y5RAfHLEHwGvoOfU3wFDy1aW9mtaya1DKLmXKlUigwuInEHYleCoLnA"
    "sGmtg9grLpSFxw6dOxnXM2AwXy6XEwOW9AMmF8BV3G/x2Cy8xx9riwbiDuIkJrwjPM/KZ5sG"
    "5t92kBebLTp7dC8GY/AWfBFiE1gU3mIGlUols8E2iu7jjLvLvldgR+AH7QM4DhXxm1gfTeWP"
    "0fRK1h18wKvox81Lw6H+aeS9FOtgfD2rwRsInEkTFnFptYTc10Ld02KxmPmINrHlE2nCIjhC"
    "F7nbvOZTvdOXZ2nyhku9xTRhEeQ9imswhg7DIKzjvG6DOxx2Pw1fG15KExVhObEBZtai5DiW"
    "HASt1OT/wab+krCDIzAwczVAs8vcYFTvDmTJtvM1wMtwihs8j17TvI8IPRgyA/TsRjTsXDdv"
    "g/AeDLbQi+PcwM7NgH0TIP4NJneTYZeXAWbUQYg/Ae/Q6GPJsMujB21/5RCEH4Ax5pbHYsIO"
    "9mYAwR6EN8CI3cfxxIDtgCj6SPMP9yT/WvOozUbDR/ACo1ubWPOac9GwYz9N06RarSrpui4p"
    "ylL0HMNi1ar+x5qu714xDkpGp7dPURRJVdUkrmkqlCXpF1KpYNTIER0eAAAAAElFTkSuQmCC")
nav3_index.append('nav3_back')
nav3_catalog['nav3_back'] = nav3_back
getnav3_backData = nav3_back.GetData
getnav3_backImage = nav3_back.GetImage
getnav3_backBitmap = nav3_back.GetBitmap
getnav3_backIcon = nav3_back.GetIcon

#----------------------------------------------------------------------
nav3_contrast = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAASFJREFUSEvNlUEOgjAQRYG1J2GnNzDs4AZuPYbn4TBwCYgegZ3B"
    "6jzTEpHSgrHGl/yElOlMOzNto38gEe1ER9FJi2/G+PcxGxHOLqL7jM4ibLBdxV5kddz3vWqa"
    "RpVlqfI8v8VxzDi2zFnEQXQVTZyjd6qqUmmaKvnHHOY6YRWzzpGNrutUlmUmyOxOyKMr30/N"
    "QRC9E3xYa0KxrE5f5aKua5UkCXb4GkG7eVePfBRFwS7orlEL09NWh+/yQXdp261ogIMzcWaT"
    "j7ZtTQB8DizKP/LBOdG2ozp8M4CxHQUInqLgRQ7eprCoDi5cBw2CXxUQ9LIzrL2ub2uuawOr"
    "sKbrGw+OgTxSLDpiEkjr4yfzFdqNnubg4AzxzdikFX9MFD0Aw8HVrGb5SeoAAAAASUVORK5C"
    "YII=")
nav3_index.append('nav3_contrast')
nav3_catalog['nav3_contrast'] = nav3_contrast
getnav3_contrastData = nav3_contrast.GetData
getnav3_contrastImage = nav3_contrast.GetImage
getnav3_contrastBitmap = nav3_contrast.GetBitmap
getnav3_contrastIcon = nav3_contrast.GetIcon

#----------------------------------------------------------------------
nav3_filesave = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAABHNCSVQICAgIfAhkiAAAAAlw"
    "SFlzAAACXwAAAl8BvoUoWgAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoA"
    "AAPeSURBVEiJtVZNa11VFF17n3OTlzRNk7ZobasUpRRttShNxCIIFhT8AB1YodiBTpyIHzMR"
    "q1J/geBAcFhURHGig86EOChKBlVLW8W2Vhqr+Whi09j3cvfey8F9X0mboKgHDufec89be6+1"
    "1+Y8IYn/c2QAWHvwjf0A7iKoAtlMoIpKEGA0z7JrxrJ3AqAAQfIHACGqR+ePHD4vA88ceofk"
    "bnrUBRx66cqXI7t2rAEJdJETMkAKQKK9TVbnUKX09bnA+3L3rABSDPTPp6J4KAMY3T2QdPvG"
    "Dbd8cnriYn//YAyuu1GqABSKACQBAcFqv8UHrIIKwCBqE9PcKpDHR3b2fTF+Ok2pvpBB9uzc"
    "csPwy/sf3vbpq+/OzD/x4kJt123m5lqaq5mJmambS+mubqYWgfCAh1erBzwM4/WzMjxxMb9y"
    "4NHasZNnbdIxnAnI5yfO/fzVyffmPOVQCMJD3F3cHWcuTqVfp+cSSSFDEY6ggBFNGR0ksWlo"
    "gCrEd5OXr+549s06c4HejeuRAchc6NVZZ6ki64UU8xALiofL2Pj3PT9+eyK3ZGlVv1WIVkVG"
    "9464CoKqNClCVBUAM0hhsACZqJJFpJm9SXjIU/v2NsoHRkorXctwCXdEVLJEGDwC4YSAmDh1"
    "QSBVsbptqmQUIBOILAxxN3F3eLh8eHSsdvL4idR2U9tYHRYAcN/9e0wFAQi6PjCDEAQLkkkU"
    "BURb+ouby77RO8t77rg1IiDBAD0QjGYNAhFBBjDYV3Dm/HRlqe5Gq6zHgqSSyH82GnrpSp1m"
    "puEuOSUZXtPPIFlJ45UsJFpyhQcW3WFBqfqtYtuSSBhRgNQorfxo7PgCSAkyAdQIauX3TlZs"
    "NllnqeQKEGC785s1IIRkBqkIw++TM7MkMhiZwYzlnP/hyAAFEQUJBZhIJgQT+e/BwWaRSWYS"
    "ClLBNvh/MjKBVh8oQQWhAJDCVSL+HgMBLBV+zWan0aLAMjm2S33w3j23a73eEGssggLknEmy"
    "ck7X/OnML3oqrZ++LoMm8DWZ9vb2yGOPPMiPP/gsfXNsPAHAweeeLrfcvDmsLGFmMDOUZYnJ"
    "S3OBhevkD1BXY27uYHRcZ+ZwM7g7zBzWfMYqt+KqAdwM0fXjpcDWfsbK+CsHEAHNDIyuAOZL"
    "gN2tYrlKBBWVYkUG7oiuxqwYdIDNDL4CA1FVCBpZVN+qrVv7dv2P+aXeJ8TMsHnrJpq5iwJ9"
    "/TV2A7ckk/a9X5U2FYX2Dq09rykdEpIYOPD6a4x4nqQ1T+Kmy1MbBntSDrK6kpvrskmSsthY"
    "xG+btk1ApPp3oXpBU3py/sjhhb8AKRQO9CWsy6YAAAAASUVORK5CYII=")
nav3_index.append('nav3_filesave')
nav3_catalog['nav3_filesave'] = nav3_filesave
getnav3_filesaveData = nav3_filesave.GetData
getnav3_filesaveImage = nav3_filesave.GetImage
getnav3_filesaveBitmap = nav3_filesave.GetBitmap
getnav3_filesaveIcon = nav3_filesave.GetIcon

#----------------------------------------------------------------------
nav3_forward = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAp9JREFUSEutlruLE0Ecx2cFvWhARUJiLFREi8vmua88LrezOZCL"
    "yWXtVBAEsbax8i+wEbSzsLCxtLLQRhA7K5tDLOwUQSx8nGhxeqLf32QmmVtHTRYDH7Lze313"
    "5rf5bRh9XLduEZ43/g4Cx2o0apbj1MTaWR7udMPYdnm8CgbAdcJR1vMawq+gPBGv8sI4y/L5"
    "POO8tSuK2guSTBS1xLXH11ZQ7B74DH4m+AYee3x0oct5dpo7rgP7ZZefPssKhQIJLPR6rYwi"
    "4KeOIvmBVuxfPPf52pLKR/GLsP2AwDm1g4kAAjtwvtWSZ2UThS+BM7j+PrZBQN+Bz4cuHB9l"
    "Qhpw1/HmdK0JtPnJAzC+nDr/BxAoFotCwOPxdXOQ4BN4Db5qthmQO0BTizAkkzcgeg3HZke9"
    "zm46wuUozKJHXdjvwL+ViDcgBZBwJeF8FvDBcdV4E3iEI8S9S+QlgEAul2NYPNIcL5p8NW8q"
    "qoMbOIHYV1qeAQjUlwYWFh+kcQtH0DQV1MGRHkHsDA8EBPBz3qcZ75sK6jR5/xDi1rWcvwAB"
    "zBhKEAb8SM6biupQDHp2U3LDwMSHmnXawX4l4PPBoqloOtqZIPAsVu8OqQcbJIDmFszBaSAB"
    "1xKPKYo/IQE075g5OA3tjO9DgIYdzusqCeB8R+bgNMgjIgGc/WEIYBrGt83BaWhnmk0IqGGH"
    "4rcg8oXeBeaEeSEB35oMOzT4IATegIeYOXvMSfMgBdQOyIj5EkKAht7dTrSy9/ekeZA90AUI"
    "iPQhQI/tOl1vT5oH2QMadmFIL3160Y/BiC5D4Cl4H/KO/BOwHbopzukF/2ef2AF9ymUblFml"
    "Uma2bbNSqcScpf4OjJG8stM3USotYl1h1WpFfCu77qvVqmLNGGO/AGh8YNRvhAXSAAAAAElF"
    "TkSuQmCC")
nav3_index.append('nav3_forward')
nav3_catalog['nav3_forward'] = nav3_forward
getnav3_forwardData = nav3_forward.GetData
getnav3_forwardImage = nav3_forward.GetImage
getnav3_forwardBitmap = nav3_forward.GetBitmap
getnav3_forwardIcon = nav3_forward.GetIcon

#----------------------------------------------------------------------
nav3_home = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAZVJREFUSEvt0s8rRFEUwPFJKdn5D8ggv2M7M7GwlZ2lEmG2ysZO"
    "WbCd0pP8KAuFspJZTVZEQ5kykgXlx0KRhfzMr+/xzsub65lmJgtqTn163XPPPfe9+64v0wiF"
    "JxuwihX4Nf07QcN23OIOj7hGi07LfAGKXYp0Kn0EByxZPIgXHKEKTTjFE7qljmc13l0Snw3S"
    "Bc0LKZzSBevBsCVNZmGhBtt4wxhqkfkGNC+hKKbFM2jEsY5FMhS26ngu6nhTn45EoG9CuxnB"
    "ZCUO8YohtOEG7gbiEgGMQr7EPee9AROtuIL80A704hnuxW736EQXHjQnvh8RyR7IjztDM+Rs"
    "zTfzIl86jBDkqySXDPRbdmPOW67YuE7soBzLOs6G/Cu5Zfs6jkhvefNpTUjTMmzpOBcxmpby"
    "jOp4TjaIYAT1OIG5KFsHqICcypJ9TgQDuTleC3KxoW2/guS5UeSQGyX33LQHr3rheYt+2iCu"
    "JSlBXo7Uq17kN/ijG1wYRY78BnaQ/P8brCHuYV5LUoK836hz7GLBrvL5PgBCeDcVgEdDxQAA"
    "AABJRU5ErkJggg==")
nav3_index.append('nav3_home')
nav3_catalog['nav3_home'] = nav3_home
getnav3_homeData = nav3_home.GetData
getnav3_homeImage = nav3_home.GetImage
getnav3_homeBitmap = nav3_home.GetBitmap
getnav3_homeIcon = nav3_home.GetIcon

#----------------------------------------------------------------------
nav3_move = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAQlJREFUSEu1kzEOwjAUQ3swDgATOzssSAyIAQYWRiZWNgaOwP2C"
    "HMnVj3FoI7XDq1In307y2y6lNCtWrLHcnhNQ/R9WdMD4+HhnWkKsqND89vpkWkKsGFHz1hAr"
    "ktP9mRabQ2a1u/TmGFPHGq2LWNEBMwZgrPM1rOh2NRRQO8mPgHt1Buv9tb8WjHUeuutJ8cKG"
    "uoAhUOMa3w9ork10uyXxVPwINCQ/ormCYi5WMOdqYkheOHsAiCGTXxFhCAqiPgbUqDkoFgEs"
    "cAFxt+5U0NUcFC9k1h+txlBADSsS7ApmgE0E8SOo7ZxYMcLG05y4hjqsqGjIWHNgRQdDWsyB"
    "FWvAuMUcWHE6UvcF5OFXiwd3+VIAAAAASUVORK5CYII=")
nav3_index.append('nav3_move')
nav3_catalog['nav3_move'] = nav3_move
getnav3_moveData = nav3_move.GetData
getnav3_moveImage = nav3_move.GetImage
getnav3_moveBitmap = nav3_move.GetBitmap
getnav3_moveIcon = nav3_move.GetIcon

#----------------------------------------------------------------------
nav3_subplots = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAABHNCSVQICAgIfAhkiAAAA+1J"
    "REFUeJzllX9I3GUcx1/P57mv+lXvXHqet7mltmW3dUxd6CxbCFYjaqwxKTYY9c+I1qDRInA2"
    "0mZQG6zBQon9EUGMYKxFi9BBNGgbK7o4ZhblLtDlneg58w517PT77Y/z9C5//DH6rzc8PPB5"
    "vt/36/k8n+cH/N90DLgIrL9Xg0/8fr/d0NBgA4tbEzaHl4hntgtAQcpQpbv7/X772rVriAgA"
    "tm1jWRbd3d2cvnCa/LZ8VjlWMd49jq/fx9GjR+f/jUajAJw/f57W1tYgUAPgSAe4XC5CoRDD"
    "w8PzEICTXSfR7Zo1s2sokRKyG7O5/ONlPB97qKurAyA7OxutNY2NjXg8nuqRkREWAVIzvn37"
    "Nv39/fNxdUcRfilMmPB8LD8/n1u5t0gkEhlrvGfPHlLmSwJmZ2dJJBKEsx5k1HwIw9Csf3kX"
    "PkNjOARjrs8yNIZDYxtCYZ5BYZ6B32tSlDWVAVwyAwARlTR0CNFQgNH+ACIKLQrRCi2CiOKJ"
    "7TtYu6WGGRu0ZJR0McCyLG7cuEFnZyfmA49R8fTD/HzuBFM3v6egoIDR0VHuc1cBMBL5Ca/X"
    "y+c/XOTm9mZeONhC66FXccrd5QGDg4N0fd2FHBfi70UwHMJv335GR0cHPp+Pnp4ecr1PAhD5"
    "3UVzczPj4+O88eZb7H29lbHcQaZ3SnKjTtMAXJU0/w8nqyZxvuNko3vj3BJpAEzTxLIsDMOg"
    "cssGKrdswDRNRATTNJmMT+DQgtPhZFP1JpzHnVBJB7B7IYPVVKsXFW7DTaku5RcVwnAk+bFY"
    "jL6+PoqLixkJfgVAeXk5vb29C0vhEBzKQYVZgdvlJt4Wh708vgCI8LbdYl+JnYkxVjaGKOYB"
    "bW1ti4qXrvb2drQWojNRZuwZxnrG4Cw3gY/Sa3C1yCwi50wOgdoAIhsoLcxZ0ThdWhRqCi6d"
    "vETsbAzgEDCZXgPcbjcHDhxgW3Qbm2rqERFWr72fYDC4rHEwGMS7rgzRwr5nD7L/kf2poUn4"
    "1y6CZEFra2v5Y9qFaEXL8S6e2bGT4b8GlwR415Xx/qdfIFp4blczrrt/Z4wvAqQkohARNm/d"
    "xrkrvyYPmChQc00UtpL53lIKewmfZQFaFFoLCQtOHHmNoYE/USp5UpUCMy+fV44co3JzDTYK"
    "y07e1alvlgSEw8nLzOPx8FRVFV6vB4BTVpSyM4X4i/2U55QTmg5x/YPrFN8ZYWthpuHExAoZ"
    "RCKRLwOBwPNNTU3E43Hi8ThKKfbt3sep9lNY71oMGUMMfTNEXVYdPp+PgYGBuaxURr+cqvLy"
    "8sZZ7rV6FFsdViu+aPX19Tbw3YqUFZR6k0vu1eA/1z/VL1jv1/scBwAAAABJRU5ErkJggg==")
nav3_index.append('nav3_subplots')
nav3_catalog['nav3_subplots'] = nav3_subplots
getnav3_subplotsData = nav3_subplots.GetData
getnav3_subplotsImage = nav3_subplots.GetImage
getnav3_subplotsBitmap = nav3_subplots.GetBitmap
getnav3_subplotsIcon = nav3_subplots.GetIcon

#----------------------------------------------------------------------
nav3_zoom_to_rect = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAABR1JREFUSEutVGtMk2cUJouJJILOG0JpQegNENzcdFMRskk0kXiB"
    "6XQqZG7BuelwOieTuIkM1IkWuSioXJQWCpaWgqWWQeUuAtLa0lY7oEBLuXhDf7Dwgx/P3q/7"
    "foww3HQ+yZOcnPd8z7m833ucXgW7o/c4c3n8cIYnU+Ds7CycOXOmcOFCN8Gy997fvXz5B850"
    "2OuBx/Pb68vm2GP27UdJaRm0uk7oOo2oUChx+MhR+AUseerhwThIh/93fH/02Cy3Re7C7bui"
    "odGb0GcbxJ2OTsjVLSiraUFjux69VjsedFnw9bffgcHwVActfXcW/fm/w5PJEh46Gk8E+pBX"
    "fhsp0lqUa7pwz/oEGttTVOosSJU3IFtWDVNXL5LPCuDpyVTTn78c3ot9IiK2bkeH/gESCxRQ"
    "mmzoGR2D5fkfsLygSWzKV9s9jJMiJdq0RsTsO4DFPr4vH9f1whszmCwvg7rxLpKECkiNNtRa"
    "n6LBPorGoedoGnrhYCMh5asj3SjMQ0gghTS1asHm8uy7oj6f/uJ9fdmrPiNzL7qpxqlqLURG"
    "O0q6RlBqeQxZ7xOU9f1FypYSH3VWaLLjfL0RObIq7I89BCaT9QktNxUeDIYgPSsHPxdU4kyT"
    "Genaflw2DiL34TDyzcO4Zh5xkLIp32XTIASafvxYY8DBnAoUFMvgwfAU03JT4e7uoRLLlIi9"
    "rkKc2oSktl6c1w0grdOOiyTRJcIMgx0p921IuNePw8092FNlRGRJKyIzyyGrvA2Wl7eelpsK"
    "NoerLyIJvshT4ZtbesTf6UZSRz/OaK34lYie0dpwithJGitJ0If41l7E1pkRJdcgIrsSpSTB"
    "IneP6RNQ/3NuoQx785SIuanFseZunCaCAv0AqZx0YBpyMJN0kqq3OxLHk5i9ivuIylGiSKYC"
    "efHTJ3BzW5R+POkcThSSDio0OEk6yCDjySUzF3U9grjnsYOUnUd81FliSzcOKLSIE6mQkpmL"
    "BQsWTn8HZC18vC48AjkSJfZL7iCp6XdcJTOXEFGF9RlUA6OoIqwkNuW7ahhAMok5IGtFtuQW"
    "tu7cQz243bTcVLDZnJB58xdMCLKFOFdajYTqTlzR9EHePYJaItwy/MJByqZ8V8nZSfIHJd+o"
    "wcX8EvLQ2I83bYn855Xh4uISxvMPGo+JjcN1SSVyxRW4IK9DesMDlBrIg+t9hHb7MwcpW0p8"
    "GY0PcU5ej7ziCqwMCQOL5bWNlpsMSpzNDxi7cKUQDXe1sJBFJlHcxsVrEhSrGlBQ24E6sx3G"
    "wWcwDY2ioWsIwjoNxFVNyC6QYs3aDZg9Z04iLTcZlDjHL3A8lYg3t+tgtQ+jtvkeTpy9BJYP"
    "B1FfHcYVURmKylSQquohrap32Hmkwy9jj4HrHzTh7b04jpabDErch+s/lnpZhKY2HWyDIzCY"
    "LRDLf8O6jdvw9ty56fPmzU9nevuOvrMiGCHrNiN0/RYs+zAUPhz+ONm64lWr1wTScpNBiXv5"
    "cMcTUrJQ09jmqNzwsAfF5dXYELEDrq6uWXSoU9GNshkrVwcvDf1obXRwSGh0wJLApWfPp82g"
    "j6fCzz8g2I3hNbbvSCLyyQW13zehXWdyVB4euZOIz87i+/m/RYe/Okj1P6wO24zjpzMhkipR"
    "qlDjErnQ9Zs+dYjTYa8P0kGgN5s/EfdLGpLTchF/KgNrwjbChYzlf1X+d/D4fnFsftDEph0x"
    "CFoR4hCnj94cyLMOcXFxzScb8CcOl/dmKndycvoTX2VtcSgrH/YAAAAASUVORK5CYII=")
nav3_index.append('nav3_zoom_to_rect')
nav3_catalog['nav3_zoom_to_rect'] = nav3_zoom_to_rect
getnav3_zoom_to_rectData = nav3_zoom_to_rect.GetData
getnav3_zoom_to_rectImage = nav3_zoom_to_rect.GetImage
getnav3_zoom_to_rectBitmap = nav3_zoom_to_rect.GetBitmap
getnav3_zoom_to_rectIcon = nav3_zoom_to_rect.GetIcon

        
        


# Test Code ------------------------------------------------------------------

class util_CreateMenuBar:
    """
    Example of the menuData function that needs to be in the program
    in which you are creating a Menu
    
        def menuData(self):
            return [("&File", (
                        ("&New",  "New Sketch file",  self.OnNew),
                        ("&Open", "Open sketch file", self.OnOpen),
                        ("&Save", "Save sketch file", self.OnSave),
                        ("", "", ""),
                        ("&Color", (
                            ("&Black",    "", self.OnColor,      wx.ITEM_RADIO),
                            ("&Red",      "", self.OnColor,      wx.ITEM_RADIO),
                            ("&Green",    "", self.OnColor,      wx.ITEM_RADIO),
                            ("&Blue",     "", self.OnColor,      wx.ITEM_RADIO),
                            ("&Other...", "", self.OnOtherColor, wx.ITEM_RADIO))),
                        ("", "", ""),
                        ("About...", "Show about window", self.OnAbout),
                        ("&Quit",    "Quit the program",  self.OnCloseWindow)))]    
    """
    def __init__(self, self2):
        menuBar = wx.MenuBar()
        for eachMenuData in self2.menuData():
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1]
            menuBar.Append(self.createMenu(self2, menuItems), menuLabel)
        self2.SetMenuBar(menuBar)

    def createMenu(self, self2, menuData):
        menu = wx.Menu()
        for eachItem in menuData:
            if len(eachItem) == 2:
                label = eachItem[0]
                subMenu = self.createMenu(self2, eachItem[1])
                menu.Append(wx.ID_ANY, label, subMenu)
            else:
                self.createMenuItem(self2, menu, *eachItem)
        return menu

    def createMenuItem(self, self2, menu, label, status, handler, kind=wx.ITEM_NORMAL):
        if not label:
            menu.AppendSeparator()
            return
        menuItem = menu.Append(-1, label, status, kind)
        self2.Bind(wx.EVT_MENU, handler, menuItem)
        
        
class DemoImagePanel(ImagePanelToolbar2):
    """Plots several lines in distinct colors."""

    # Activate event messages
    _EVENT_DEBUG = True
    
    def __init__( self, parent, tab, statusbar, **kwargs ):
        # statusbar has to be here for NavigationToolbar3Wx to discover on init()
        self.statusbar = statusbar
        # initiate plotter
        sizer = ImagePanelToolbar2.__init__( self, 
                                             parent, 
                                             vertOn=True, 
                                             horizOn=True, 
                                             lcolor='gold',
                                             lw=0.5,
                                             **kwargs )  
        
        self.tab    = tab
        self.top    = wx.GetApp().GetTopWindow()
        self.parent = parent
        self.count  = 0
        
        self.statusbar = statusbar

    def on_motion(self, xloc, yloc, iplot):
        value = self.data[iplot][0]['data'][int(round(xloc))][int(round(yloc))]
        self.top.statusbar.SetStatusText( " Value = %s" % (str(value), ), 0)
        self.top.statusbar.SetStatusText( " X,Y = %i,%i" % (int(round(xloc)),int(round(yloc))) , 1)
        self.top.statusbar.SetStatusText( " " , 2)
        self.top.statusbar.SetStatusText( " " , 3)

    def on_panzoom_motion(self, xloc, yloc, iplot):
        axes = self.axes[iplot]
        xmin,xmax = axes.get_xlim()
        ymax,ymin = axes.get_ylim()         # max/min flipped here because of y orient top/bottom
        xdelt, ydelt = xmax-xmin, ymax-ymin
        
        self.top.statusbar.SetStatusText(( " X-range = %.1f to %.1f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText(( " Y-range = %.1f to %.1f" % (ymin, ymax)), 1)
        self.top.statusbar.SetStatusText(( " delta X,Y = %.1f,%.1f " % (xdelt,ydelt )), 2)
        self.top.statusbar.SetStatusText( " " , 3)

    def on_level_motion(self, xloc, yloc, iplot):
        self.top.statusbar.SetStatusText( " " , 0)
        self.top.statusbar.SetStatusText(( " Width = %i " % (self.width[iplot],)), 1)
        self.top.statusbar.SetStatusText(( " Level = %i " % (self.level[iplot],)), 2)
        self.top.statusbar.SetStatusText( " " , 3)



class MyFrame(wx.Frame):
    def __init__(self, title="New Title Please", size=(350,200)):
 
        wx.Frame.__init__(self, None, title=title, pos=(150,150), size=size)
        self.Bind(wx.EVT_CLOSE, self.on_close)
 
        util_CreateMenuBar(self)

        self.statusbar = self.CreateStatusBar(4, 0)
 
        self.size_small  = 64
        self.size_medium = 128
        self.size_large  = 256
 
        data1 = { 'data'  : self.dist(self.size_medium),
                  'alpha' : 1.0
                }
 
        data2 = { 'data'  : 100-self.dist(self.size_medium),
                  'alpha' : 0.5,
                  'cmap'  : cm.hsv,
                }
         
        data = [[data1], [data2]]
         
        # data1 = {'data': self.dist(self.size_medium),
        #          'alpha': 1.0
        #          }
        #
        # data = [[data1],]
         
        self.nb = wx.Notebook(self, -1, style=wx.BK_BOTTOM)
         
        panel1 = wx.Panel(self.nb, -1)
         
        self.view = DemoImagePanel(panel1, self, self.statusbar, naxes=2, data=data)
#        self.view = DemoImagePanel(panel1, self, self.statusbar, naxes=1, data=data)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        panel1.SetSizer(sizer)
        self.view.Fit()    
     
        self.nb.AddPage(panel1, "One")

    def menuData(self):
        return [("&File", (
                    ("", "", ""),
                    ("&Quit",    "Quit the program",  self.on_close))),
                 ("Tests", (
                     ("Set Small Images - keep norm",  "", self.on_small_images_keep_norm),
                     ("Set Medium Images - keep norm", "", self.on_medium_images_keep_norm),
                     ("Set Large Images - keep norm",  "", self.on_large_images_keep_norm),
                     ("", "", ""),
                     ("Set Small Images",  "", self.on_small_images),
                     ("Set Medium Images", "", self.on_medium_images),
                     ("Set Large Images",  "", self.on_large_images),
                     ("", "", ""),
                     ("Placeholder",    "non-event",  self.on_placeholder)))]       

    def on_close(self, event):
        dlg = wx.MessageDialog(self, 
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()

    def on_placeholder(self, event):
        print("Event handler for on_placeholder - not implemented")


    def on_small_images_keep_norm(self, event):
        self.on_small_images(event, keep_norm=True)

    def on_small_images(self, event, keep_norm=False):

        data1 = { 'data'  : self.dist(self.size_small),
                  'alpha' : 1.0
                }

        data2 = { 'data'  : 100-self.dist(self.size_small),
                  'alpha' : 0.5,
                  'cmap'  : cm.hsv,
                }
        
        data = [[data1], [data2]]

        self.view.set_data(data, keep_norm=keep_norm)
        self.view.update(no_draw=True, keep_norm=keep_norm)
        self.view.canvas.draw()

    def on_medium_images_keep_norm(self, event):
        self.on_medium_images(event, keep_norm=True)

    def on_medium_images(self, event, keep_norm=False):

        data1 = { 'data'  : self.dist(self.size_medium),
                  'alpha' : 1.0
                }

        data2 = { 'data'  : 100-self.dist(self.size_medium),
                  'alpha' : 0.5,
                  'cmap'  : cm.hsv,
                }
        
        data = [[data1], [data2]]

        self.view.set_data(data, keep_norm=keep_norm)
        self.view.update(no_draw=True, keep_norm=keep_norm)
        self.view.canvas.draw()

    def on_large_images_keep_norm(self, event):
        self.on_large_images(event, keep_norm=True)

    def on_large_images(self, event, keep_norm=False):

        data1 = { 'data'  : self.dist(self.size_large),
                  'alpha' : 1.0
                }

        data2 = { 'data'  : 100-self.dist(self.size_large),
                  'alpha' : 0.5,
                  'cmap'  : cm.hsv,
                }
        
        data = [[data1], [data2]]

        self.view.set_data(data, keep_norm=keep_norm)
        self.view.update(no_draw=True, keep_norm=keep_norm)
        self.view.canvas.draw()        
        

    def dist(self, n, m=None):  
        """
        Return a rectangular array in which each pixel = euclidian
        distance from the origin.

        """
        n1 = n
        m1 = m if m else n

        x = np.arange(n1)
        x = np.array([val**2 if val < (n1-val) else (n1-val)**2 for val in x ])
        a = np.ndarray((n1,m1),float)   # Make array

        for i in range(int((m1/2)+1)):  # Row loop
            y = np.sqrt(x + i**2.0)     # Euclidian distance
            a[i,:] = y                  # Insert the row
            if i != 0: a[m1-i,:] = y    # Symmetrical
        return a


#----------------------------------------------------------------
#----------------------------------------------------------------
# Saved code

# class MyNavigationToolbar(NavigationToolbar2WxAgg):
#     """
#     Extend the default wx toolbar with your own event handlers
#     """
#     ON_CUSTOM = wx.NewIdRef()
#     def __init__(self, canvas, cankill):
#          
#         # create the default toolbar
#         NavigationToolbar2WxAgg.__init__(self, canvas)
#  
#         # remove the unwanted button
#         POSITION_OF_CONFIGURE_SUBPLOTS_BTN = 6
#         self.DeleteToolByPos(POSITION_OF_CONFIGURE_SUBPLOTS_BTN) 
#  
#         # for simplicity I'm going to reuse a bitmap from wx, you'll
#         # probably want to add your own.
#         self.AddSimpleTool(self.ON_CUSTOM, _load_bitmap('stock_left.xpm'), 'Click me', 'Activate custom contol')
#         wx.EVT_TOOL(self, self.ON_CUSTOM, self._on_custom)
#  
#     def _on_custom(self, evt):
#         # add some text to the axes in a random location in axes (0,1)
#         # coords) with a random color
#  
#         # get the axes
#         ax = self.canvas.figure.axes[0]
#  
#         # generate a random location can color
#         x,y = tuple(rand(2))
#         rgb = tuple(rand(3))
#  
#         # add the text and draw
#         ax.text(x, y, 'You clicked me',
#                 transform=ax.transAxes,
#                 color=rgb)
#         self.canvas.draw()
#         evt.Skip()



#     def calc_lut_value(self, data, width, level):
#         """Apply Look-Up Table values to data for given width/level values."""
# 
#         conditions = [data <= (level-0.5-(width-1)/2), data > (level-0.5+(width-1)/2)]
#         functions  = [0, 255, lambda data: ((data - (level - 0.5))/(width-1) + 0.5)*(255-0)]   # 3rd function is default
#         lutvalue = np.piecewise(data, conditions, functions)
# 
#         # Convert the resultant array to an unsigned 8-bit array to create
#         # an 8-bit grayscale LUT since the range is only from 0 to 255
#         return np.array(lutvalue, dtype=np.uint8)





if __name__ == '__main__':

    app   = wx.App( False )
    frame = MyFrame( title='WxPython and Matplotlib', size=(600,600) )
    frame.Show()
    app.MainLoop()

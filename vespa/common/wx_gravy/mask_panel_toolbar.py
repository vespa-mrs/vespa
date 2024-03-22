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
"""
# Python modules

import os


# Used to guarantee to use at least Wx2.8
#import wxversion
#wxversion.ensureMinimal('2.8')

# third party modules
import wx
import numpy as np
import matplotlib
import matplotlib.cm as cm

matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backend_bases          import NavigationToolbar2, cursors

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


class MaskPanelToolbar2(wx.Panel):
    """
    The ImagePanel has a Figure and a Canvas and 'n' Axes. The user defines
    the number of axes on Init and this number cannot be changed thereafter.
    However, the user can change the number of axes displayed in the Figure.
    
    Axes are specified on Init because the zoom and widlev actions 
    need an axes to attach to initialize properly.  
    
    on_size events simply set a flag, and the actual resizing of the figure is 
    triggered by an Idle event.
    
    """

    # Set _EVENT_DEBUG to True to activate printing of   
    # messages to stdout during events.
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
                               **kwargs):
        
        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        self.figure = Figure(figsize=(5,4), dpi=100)

        self.imageid    = [None    for i in range(naxes)]
        self.cmap       = [cm.gray for i in range(naxes)]
        self.width      = [WIDMAX  for i in range(naxes)]
        self.level      = [LEVSTR  for i in range(naxes)]
        self.vmax       = [WIDMAX  for i in range(naxes)]
        self.vmin       = [LEVSTR  for i in range(naxes)]
        self.vmax_orig  = [WIDMAX  for i in range(naxes)]
        self.vmin_orig  = [LEVSTR  for i in range(naxes)]
        self.mask       = [None    for i in range(naxes)]
        self.mask_alpha = [1.0     for i in range(naxes)]
        self.mask_cmap  = [cm.Reds for i in range(naxes)]
        
        self.paint_mode = 'add_voxels'
        self.paint_radius = 3
        
        # here we create the required naxes, add them to the figure, but we
        # also keep a permanent reference to each axes so they can be added
        # or removed from the figure as the user requests 1-N axes be displayed
        self.naxes = naxes   
        self.axes  = []
        self.axes.append(self.figure.add_subplot(naxes,1,1))
        if naxes > 1:
            iaxes = np.arange(naxes-1,dtype=np.uint16)+1
            for i in iaxes:
                self.axes.append(self.figure.add_subplot(naxes,1,i+1, sharex=self.axes[0], sharey=self.axes[0]))
        self.all_axes = list(self.axes)
 
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.figure.subplots_adjust(left=0.0,right=1.0,
                                    bottom=0.0,top=1.0,
                                    wspace=0.0,hspace=0.0)

        # for images we don't show x or y axis values
        for axes in self.axes:
            axes.set_adjustable('box-forced')
            axes.xaxis.set_visible(False)
            axes.yaxis.set_visible(False)

        self.set_color( color )
        self.set_axes_color()

        if not data or len(data) != naxes:
            data = self._default_data()
        self.set_data(data)
        self.update(no_draw=True)   # we don't have a canvas yet

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        # Capture the paint message
        wx.EVT_PAINT(self, self.on_paint)

        self.toolbar = NavigationToolbar3Wx(self.canvas, 
                                            self, 
                                            vertOn=vertOn, 
                                            horizOn=horizOn,
                                            lcolor=lcolor,
                                            lw=lw)
        self.toolbar.Realize()
        if wx.Platform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSize()
            fw, fh = self.canvas.GetSize()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wx.Size(fw, th))
            #self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
            self.sizer.Add(self.toolbar, 0, wx.ALIGN_CENTER | wx.EXPAND)

        # update the axes menu on the toolbar
        self.toolbar.update()
        
        self.SetSizer(self.sizer)
        self.Fit()

        # Always link scroll events, but only handle them if on_scroll()
        # is overloaded by the user.
        self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)


    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================

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
        m1 = m if m else n

        x = np.arange(n1)
        x = np.array([val**2 if val < (n1-val) else (n1-val)**2 for val in x ])

        a = np.ndarray((n1,m1),float)   # Make array

        for i in range(int((m1/2)+1)):  # Row loop
            y = np.sqrt(x + i**2.0)     # Euclidian distance
            a[i,:] = y                  # Insert the row
            if i != 0:
                a[m1-i,:] = y           # Symmetrical

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

    def on_paint_press(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_paint_press,     xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

    def on_paint_release(self, xloc, yloc):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_paint_release,   xloc='+str(xloc)+'  yloc='+str(yloc))

    def on_paint_motion(self, xloc, yloc, iplot):
        """ placeholder, overload for user defined event handling """
        self._dprint('debug::on_paint_motion,    xloc='+str(xloc)+'  yloc='+str(yloc)+'  Index = '+str(iplot))

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
                else:
                    # Only data in this item, so add all default values 
                    dat = { 'data'      : dat,
                            'alpha'     : 1.0,
                            'cmap'      : self.cmap[i],
                            'vmax'      : dat.max(),
                            'vmin'      : dat.min(),
                            'keep_norm' : keep_norm
                          }

                if self.mask[i] is None:
                    self.mask[i] = np.zeros(dat['data'].shape, dtype=np.byte)
                else:
                    if self.mask[i].shape != dat['data'].shape:
                        self.mask[i] = np.zeros(dat['data'].shape, dtype=np.byte)
                        
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
            
            axes.images = []
            
            ddict  = self.data[i][0]
            data   = ddict['data'].copy() 
            alpha  = ddict['alpha']
            cmap   = ddict['cmap']
            vmax   = self.vmax[i]
            vmin   = self.vmin[i]
            
            mask   = np.ma.masked_where(self.mask[i] <= 0, self.mask[i])
            cmap2  = self.mask_cmap[i]
            alpha2 = self.mask_alpha[i]
            
            xmin, xwid, ymin, ywid = 0, data.shape[1], 0, data.shape[0]
                        
            # Explicit set origin here so the user rcParams value is not used. 
            # This keeps us consistent across users.
            self.imageid[i] = axes.imshow(data, cmap=cmap, 
                                                alpha=alpha, 
                                                vmax=vmax, 
                                                vmin=vmin, 
                                                aspect='equal', 
                                                origin='upper') 

            self.imageid[i] = axes.imshow(mask, cmap=cmap2, 
                                                alpha=alpha2, 
                                                vmax=1.0, 
                                                vmin=0.0, 
                                                aspect='equal', 
                                                origin='upper',
                                                interpolation='nearest') 

            if xold != xwid or yold!=ywid or force_bounds: 
                xmin -= 0.5     # this centers the image range to voxel centers
                ymin -= 0.5   
                # Set new bounds for dataLims to x,y extent of data in the 
                # new image. On reset zoom this is how far we reset the limits.
                axes.ignore_existing_data_limits = True
                axes.update_datalim([[xmin,ymin],[xmin+xwid,ymin+ywid]])


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
    #   text,           # the text of the button (often not visible to users)
    #   tooltip_text,   # the tooltip shown on hover (where possible)
    #   image_file,     # name of the image for the button (without the extension)
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
        ('Paint', 'Set mask voxels On/Off RightMouse/LeftMouse', 'nav3_paintbrush_32x32', 'paint'),
        ('Cursors', 'Track mouse movement with crossed cursors', 'nav3_crosshair', 'cursors'),
        (None, None, None, None),
        ('Subplots', 'Configure subplots', 'nav3_subplots', 'configure_subplots'),
        ('Save', 'Save the figure', 'nav3_filesave', 'save_figure'),
      )
    
    
    def __init__(self, canvas, parent, vertOn=False, horizOn=False, lcolor='gold', lw=0.5):
        
        wx.ToolBar.__init__(self, canvas.GetParent(), -1)
        NavigationToolbar2.__init__(self, canvas)

        self.parent  = parent
        self.statbar = self._parent.statusbar
        self._idle   = True
        self._idRelease = self.canvas.mpl_connect('button_release_event', self.release)
        # turn off crosshair cursors when mouse outside canvas
        self._idAxLeave  = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self._idFigLeave = self.canvas.mpl_connect('figure_leave_event', self.leave)

        axes = self.parent.axes
        
        # set up control params for cursor crosshair functionality
        self._cursors = False
        self.vertOn  = vertOn
        self.horizOn = horizOn
        self.vlines = []
        self.hlines = []
        if vertOn:  self.vlines = [ax.axvline(0, visible=False, color=lcolor, lw=lw) for ax in axes]
        if horizOn: self.hlines = [ax.axhline(0, visible=False, color=lcolor, lw=lw) for ax in axes]

        
    def _init_toolbar(self):
        self._parent = self.canvas.GetParent()

        self.wx_ids = {}
        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.AddSeparator()
                continue
            
            bmp = icon_catalog[image_file].GetBitmap()
            
            self.wx_ids[text] = wx.NewIdRef()
            if text in ['Pan', 'Level', 'Cursors', 'Paint']:
                self.AddCheckTool(self.wx_ids[text], bmp, shortHelp=text, longHelp=tooltip_text)
            elif text in ['Zoom', 'Subplots']:
                pass    # don't want this in my toolbar
            else:
                self.AddSimpleTool(self.wx_ids[text],  bmp, text, tooltip_text)
            
            self.Bind(wx.EVT_TOOL, getattr(self, callback), id=self.wx_ids[text])

        self.Realize()
    

    def pan(self, *args):
        """ 
        Toggle the pan/zoom tool on/off  
        
        In this state, the left mouse button will drag the image around, the
        right button will zoom the image perspective in/out (not a rubber-band
        zoom).
        
        """
        # Turn associated buttons on/off as required 
        if 'Zoom' in list(self.wx_ids.keys()):
            self.ToggleTool(self.wx_ids['Zoom'], False)
        if 'Level' in list(self.wx_ids.keys()):
            self.ToggleTool(self.wx_ids['Level'], False)
        if 'Paint' in list(self.wx_ids.keys()):
            self.ToggleTool(self.wx_ids['Paint'], False)

        NavigationToolbar2.pan(self, *args)
        

    def level(self, *args):
        """
        Toggle the width/level tool on/off

        Change values with right button, up/down changes level value, and the 
        left/right motion changes the width value
        
        """
        if self._active == 'LEVEL':
            if 'Level' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Level'], False)
        else: 
            if 'Level' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Level'], True)
            if 'Pan' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Pan'], False)
            if 'Zoom' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Zoom'], False)      
            if 'Paint' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Paint'], False)
        
        # set pointer icon and button press funcs to the appropriate callbacks
    
        if self._active == 'LEVEL':
            self._active = None
        else:
            self._active = 'LEVEL'

        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''
        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''
    
        if self._active:
            self._idPress   = self.canvas.mpl_connect( 'button_press_event',   self.press_level)
            self._idRelease = self.canvas.mpl_connect( 'button_release_event', self.release_level)
            self.mode = 'width/level'
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)
    
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)
    
        self.set_message(self.mode)
 

    def paint(self, *args):
        """
        Toggle the mask paint tool on/off
        
        Turn on voxels with the left button and right button turns off voxels
        
        So, quick overview ...
        
        paint() - toggles on the specific Mode (ie. paint or no paint)
        
        From here on, we've set the button_press_event and button_release_events
        to deal with Paint mode stuff. The press_paint() and release_paint() 
        methods cope with what we do with our mouse in the canvas.  The 
        mouse_move() method calls the on_motion() method in which we deal with
        whatever mode is active.
        
        """
        if self._active == 'PAINT':
            if 'Paint' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Paint'], False)
        else: 
            if 'Paint' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Paint'], True)
            if 'Level' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Level'], False)
            if 'Pan' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Pan'], False)
            if 'Zoom' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Zoom'], False)      
        
        # set the pointer icon and button press funcs to the
        # appropriate callbacks
    
        if self._active == 'PAINT':
            self._active = None
        else:
            self._active = 'PAINT'
        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''
    
        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''
    
        if self._active:
            self._idPress   = self.canvas.mpl_connect( 'button_press_event',   self.press_paint)
            self._idRelease = self.canvas.mpl_connect( 'button_release_event', self.release_paint)
            self.mode = 'paint'
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)
    
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)
    
        self.set_message(self.mode)
 

    def cursors(self, *args):
        """
        Toggle the crosshair cursors tool 
        
        This will show vertical and horizontal lines that track the mouse 
        motion, or not. Only displayed in width/level and paint modes.
        
        """
        if self._cursors:
            self._cursors = False
            if 'Cursors' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Cursors'], False)
        else: 
            self._cursors = True
            if 'Cursors' in list(self.wx_ids.keys()):
                self.ToggleTool(self.wx_ids['Cursors'], True)


    # def press_pan() is in base class and ends with self.press(event) 

    def drag_pan(self, event):
        """
        The native drag_pan method ends with dynamic_update()
        We overwrite it here, call the original, then call mouse_move() method
        
        """
        NavigationToolbar2.drag_pan(self, event)
        self.mouse_move(event)

    # def release_pan() is in base class and ends with self.release(event)
        
 
    def press_level(self, event):
        """the press mouse button in width/level mode callback"""
        
        if event.button == 3:
            self._button_pressed = 3
        else:
            self._button_pressed = None
            return
        
        x, y = event.x, event.y
        
        # push the current view to define home if stack is empty
        if self._views.empty():
            self.push_current()
        
        self._xypress = []
        for i, a in enumerate(self.canvas.figure.get_axes()):
            if (x is not None and y is not None and a.in_axes(event) and a.get_navigate() and a.can_pan()):
                self._xypress.append([a, i, event.x, event.y])
                self.canvas.mpl_disconnect(self._idDrag)
                self._idDrag = self.canvas.mpl_connect('motion_notify_event', self.drag_level)
        
        self.press(event)
 
 
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
            vmin  = vmid - (vwid/2.0)
            vmax  = vmid + (vwid/2.0)
            
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
        self.release(event)
        self.draw()

 
    def press_paint(self, event):
        """the press mouse button in paint mode callback"""
        
        if event.button == 1:
            self._button_pressed = 1
        elif event.button == 3:
            self._button_pressed = 3
        else:
            self._button_pressed = None
            return
        
        x, y = event.x, event.y
        
        self._xypress = []
        for i, a in enumerate(self.canvas.figure.get_axes()):
            if (x is not None and y is not None and a.in_axes(event) and a.get_navigate() and a.can_pan()):
                self._xypress.append([a, i, event.x, event.y])
                self.canvas.mpl_disconnect(self._idDrag)
                self._idDrag = self.canvas.mpl_connect('motion_notify_event', self.drag_paint)

        self.dynamic_update()
        self.press(event)
 
 
    def drag_paint(self, event):
        """the drag callback in width/level mode"""
        
        # prep new 'last' location for next event call
        self._xypress[0][2] = event.x
        self._xypress[0][3] = event.y

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

        if iplot is not None:
            if self._cursors and self.vertOn:
                for line in self.vlines:
                    line.set_xdata((xloc, xloc))
                    line.set_visible(True)
            if self._cursors and self.horizOn:
                for line in self.hlines:
                    line.set_ydata((yloc, yloc))
                    line.set_visible(True)    
        if (self._button_pressed == 1 or self._button_pressed == 3):
            self.toggle_mask(iplot, int(event.ydata), int(event.xdata), width=self.parent.paint_radius)
            self.parent.update_images(index=iplot)
            
        self.mouse_move(event)
        self.dynamic_update()

 
    def release_paint(self, event):
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
        self.release(event)
        self.draw()

 
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
        cursor2 = wx.StockCursor(cursor)
        self.canvas.SetCursor( cursor2 )


    def press(self, event):
        
        xloc, yloc = self.get_bounded_xyloc(event)        
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
        elif self.mode == 'paint':
            for line in self.vlines + self.hlines:
                if self._cursors:
                    line.set_visible(True)
                else:
                    line.set_visible(False)

            if self._button_pressed == 1:
                self.parent.paint_mode = "add_voxels"
            else: 
                self.parent.paint_mode = "erase_voxels"
            
            self.toggle_mask(iplot, int(event.ydata), int(event.xdata), width=self.parent.paint_radius)
            self.parent.update_images(index=iplot)
            self.dynamic_update()
            self.parent.on_paint_press(xloc, yloc, iplot)
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
    
    
    def release(self, event):
        # find out what plot we released in
        iplot = None
        for i,axes in enumerate(self.parent.axes):
            if axes == event.inaxes:
                iplot = i
        
        # get bounded location and call user event
        xloc, yloc = self.get_bounded_xyloc(event)
        if   self.mode == 'pan/zoom':
            self.parent.on_panzoom_release(xloc, yloc)
        elif self.mode == 'zoom rect':
            pass
        elif self.mode == 'width/level':
            self.parent.on_level_release(xloc, yloc)
        elif self.mode == 'paint':
            self.parent.on_paint_release(xloc, yloc)
        elif self.mode == '':
            # no toggle buttons on
            self.parent.on_select(xloc, yloc, iplot)
        else:
            # catch all
            self.parent.on_select(xloc, yloc)


    def leave(self, event):
        """ Turn off the cursors as we move mouse outside axes or figure """

        for line in self.vlines+self.hlines:
            line.set_visible(False)
        self.dynamic_update()            
            

    def mouse_move(self, event):
        # call base event
        NavigationToolbar2.mouse_move(self, event)

        # The following limit when motion events are triggered 
        
        if not event.inaxes and self.mode == '': 
            return
        if not event.inaxes and not self._button_pressed:
            # When we are in pan or zoom or level and the mouse strays outside
            # the canvas, we still want to have events even though the xyloc
            # can not be properly calculated. We are reporting other things 
            # that do not need xylocs for like width/level values. 
            return
        
        # find out which axis we are in and what x,y location
        
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
        
        if self.mode == 'pan/zoom':
            if event.inaxes and self._lastCursor != wx.CURSOR_HAND:
                self.set_cursor(wx.CURSOR_HAND)
                self._lastCursor = wx.CURSOR_HAND
            if (self._button_pressed == 1 or self._button_pressed == 3):
                if iplot is not None:
                    self.parent.on_panzoom_motion(xloc, yloc, iplot)
        
        elif self.mode == 'zoom rect':
            if event.inaxes and self._lastCursor != wx.CURSOR_ARROW:
                self.set_cursor(wx.CURSOR_ARROW)
                self._lastCursor = wx.CURSOR_ARROW
            if (self._button_pressed == 1 or self._button_pressed == 3):
                if iplot is None:
                    pass
        
        elif self.mode == 'width/level':
            if event.inaxes and self._lastCursor != wx.CURSOR_BULLSEYE:
                self.set_cursor(wx.CURSOR_BULLSEYE)
                self._lastCursor = wx.CURSOR_BULLSEYE
            if self._button_pressed == 3:
                if iplot is not None:
                    self.parent.on_level_motion(xloc, yloc, iplot)
        
        elif self.mode == 'paint':
            if event.inaxes and self._lastCursor != wx.CURSOR_CROSS:
                self.set_cursor(wx.CURSOR_CROSS)
                self._lastCursor = wx.CURSOR_CROSS
            if iplot is not None:
                self.parent.on_paint_motion(xloc, yloc, iplot)

        elif self.mode == '':
            if event.inaxes and self._lastCursor != wx.CURSOR_ARROW:
                self.set_cursor(wx.CURSOR_ARROW)
                self._lastCursor = wx.CURSOR_ARROW
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
            
  

    def toggle_mask(self, iplot, yctr, xctr, width=1):
        
        xmax, ymax = self.parent.mask[iplot].shape
        
        val = 1 if self._button_pressed == 1 else 0  
        
        xtr = 1 - (width % 2)
        wid = int(width / 2)
        
        ystr = yctr - wid + xtr
        yend = yctr + wid + 1
        xstr = xctr - wid + xtr
        xend = xctr + wid + 1
        
        xstr = xstr if xstr > 0 else 0
        ystr = ystr if ystr > 0 else 0
        xend = xend if xend < xmax else xmax
        yend = yend if yend < ymax else ymax 
        
        self.parent.mask[iplot][ystr:yend, xstr:xend] = val
        

    def get_bounded_xyloc(self, event):
        if not event.inaxes:
            return 0,0

        xloc, yloc = event.xdata, event.ydata
        
        # now we bound these values to be inside the size of the image. This
        # is needed, because a pan event could yield negative locations.
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
        xmin,xmax = x0+0.5, x0+x1
        ymin,ymax = y0+0.5, y0+y1     # position swap due to 0,0 location 'upper'
        xloc = max(0, min(xmax, xloc))
        yloc = max(0, min(ymax, yloc))

#         # if the data was padded to create a square image for clarity in display,
#         # then we need to adjust the x,y display values to data array values.
#         for i,axes in enumerate(self.parent.axes):
#             if axes == event.inaxes:
#                 iplot = i
        
        return xloc,yloc
    
    
    def dynamic_update(self):
        d = self._idle
        self._idle = False
        if d:
            self.canvas.draw()
            self._idle = True


    def draw_rubberband(self, event, x0, y0, x1, y1):
        """ 
        See backend_wx.py  NavigationToolbar2Wx() method
        
        We are not using the zoom() button in this incarnation, so
        we do not need the rubberband functionality, I think. 
        
        """
        pass
        
    
    def set_status_bar(self, statbar):
        self.statbar = statbar
    
    
    def set_message(self, s, index=0):
       if self.statbar is not None:
           self.statbar.SetStatusText(s,index)
    
    
    def set_history_buttons(self):
        can_backward = (self._views._pos > 0)
        can_forward = (self._views._pos < len(self._views._elements) - 1)
        self.EnableTool(self.wx_ids['Back'], can_backward)
        self.EnableTool(self.wx_ids['Forward'], can_forward)


    def get_canvas(self, frame, fig):
        """ was in NavigationToolbar2WxAgg, needed for configure subplots """
        return FigureCanvasWxAgg(frame, -1, fig)
    


# ***************** Embed Icon Catalog Starts Here *******************

icon_catalog = {}
icon_index = []


#----------------------------------------------------------------------
nav3_crosshair = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAIBJREFUWEftlVEKgCAQBfcgXau/TtB9OkF/XbJ8sBsimlToFryB"
    "+RJxQGXl70yqG5vqBgMYwAAGrGpXhuBYEGvNwUF7Qaw1xwJSXgVgotmDuhL3XQtYgum+nHPw"
    "xD3gDrWA5lhAzi4B7t8wBvcN3bAH5QYDGMAABmCiPZ5qH0DkAEnpVB1zVLSlAAAAAElFTkSu"
    "QmCC")
icon_index.append('nav3_crosshair')
icon_catalog['nav3_crosshair'] = nav3_crosshair
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
icon_index.append('nav3_back')
icon_catalog['nav3_back'] = nav3_back
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
icon_index.append('nav3_contrast')
icon_catalog['nav3_contrast'] = nav3_contrast
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
icon_index.append('nav3_filesave')
icon_catalog['nav3_filesave'] = nav3_filesave
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
icon_index.append('nav3_forward')
icon_catalog['nav3_forward'] = nav3_forward
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
icon_index.append('nav3_home')
icon_catalog['nav3_home'] = nav3_home
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
icon_index.append('nav3_move')
icon_catalog['nav3_move'] = nav3_move
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
icon_index.append('nav3_subplots')
icon_catalog['nav3_subplots'] = nav3_subplots
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
icon_index.append('nav3_zoom_to_rect')
icon_catalog['nav3_zoom_to_rect'] = nav3_zoom_to_rect
getnav3_zoom_to_rectData = nav3_zoom_to_rect.GetData
getnav3_zoom_to_rectImage = nav3_zoom_to_rect.GetImage
getnav3_zoom_to_rectBitmap = nav3_zoom_to_rect.GetBitmap
getnav3_zoom_to_rectIcon = nav3_zoom_to_rect.GetIcon

#----------------------------------------------------------------------
nav3_paintbrush_24x24 = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAZZJREFUSEvtlD1Lw1AUhpvP5rNphRgwpQ1W12qkmyk4SQepkFn/"
    "g4L9A26Ko5PgFtw7dVdwcHQRQReHDuKgCAW/qO89pv0B3jv2gct9uW95z83JSQszGM1mUymV"
    "Stu6rp8Vi8VL0zSzNE3V3ObDdd01VVXvIMdsSZL0iWIxmbyEYTgvy/ILJIWzhac4Zp4QDMPY"
    "xzYNx+1H9XrdZx43CAtQ4AhyWgD9P2UeN3iJG2jNE/rfxv6G9YVi53EcG/lP/o9t2yssFJLd"
    "+gq3bkVR5JLJS6VSWVQUZQhJLUGhVxRYJ5MXhPsIv4ecvNAhxnGVTF5833cw6zeQk5s/lsvl"
    "JTJ5aTQaGsIHkBSOp7hF+AKZvHQ6HRk9ziApHIWugyCYI1MElmWdYJuED2q1mk2GCDCOB9go"
    "XNO0iyRJdDJEgOnYxfaDRV9nt9tVyBABbr6JEfyAHKNFh71eT/pzBIDwFkbwHfIb4XtZxt6v"
    "IDzPW0b4M/s/R6Gd/FgMmHUD8/2AAiPHcbbyY3FUq1UXk9LHB9TOj2bkFAq/C4FSOxKvVq8A"
    "AAAASUVORK5CYII=")
icon_index.append('nav3_paintbrush_24x24')
icon_catalog['nav3_paintbrush_24x24'] = nav3_paintbrush_24x24
getnav3_paintbrush_24x24Data = nav3_paintbrush_24x24.GetData
getnav3_paintbrush_24x24Image = nav3_paintbrush_24x24.GetImage
getnav3_paintbrush_24x24Bitmap = nav3_paintbrush_24x24.GetBitmap
getnav3_paintbrush_24x24Icon = nav3_paintbrush_24x24.GetIcon

#----------------------------------------------------------------------
nav3_paintbrush_32x32 = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAAj1JREFUWEftlT2LE1EUhieZfGfyMfk0JI0xjRAjBiQkFjYKERRs"
    "VFAwBAQbl1TpshArU4ZNJakSJJXsNuofWPQHaGFjoSDbqIgoirJEn3tNQqL1mWofuJx35gy8"
    "M3fuOcc4wimm06nR7XYtVng2my3uClOpVHzxePyWz+fb9Xq9B6Zpzj0ez9dMJnNq8YgcyWTy"
    "EsZvkb/XVzQavUGUJRQKXXO5XIfIDfNgMLjb6/WQglSrVYtt/oDcMHe73b9s2y6iZYnFYlcI"
    "G+Zq+f1++ZOXTqfdvMBd5IY5v2OeSqXOouXgcNUJj3O53BninLUypxK2J5MJl0JYlnWSEvuI"
    "VAfvArtwh3/+PBAI7FENFweDgX5OBGo6j/k75HLLv7O6kUjEIsrCl9o0mFfI1f9Wi3tPEolE"
    "AC1Ho9EI0Wj2kRvmnPZHtVrNhxbFxHyP+K/5TqvVcqPl6Pf7LjraQ+TKmJOuzLfrdVUIwtBm"
    "7xPWv/yQA7dFlAdz1WTWa/wn5jdVThzK7Tq1vRow6G/U+GWdlAaj8xj+QC7NP9H5zumkNIVC"
    "4TRb/RmpzZl0B9lstqKT0hSLxeN0uffIpfkbduOETkqTz+dTGL5GanO620vm+TGdlKZcLocx"
    "f4HU5tT4fqlUsnXSATxMsafE5Zc/o+3KDxZFp9NRXU4Nbm3Oi0ybzaZXJ6UZDodGOBx+gNTm"
    "vMhOu902ddIJ6GhXCbqvo/uj0ehvwikwvc3B+8KMvzcejxd3HUSZcthkR+kR/2EYfwBbDZ2+"
    "f2U0zgAAAABJRU5ErkJggg==")
icon_index.append('nav3_paintbrush_32x32')
icon_catalog['nav3_paintbrush_32x32'] = nav3_paintbrush_32x32
getnav3_paintbrush_32x32Data = nav3_paintbrush_32x32.GetData
getnav3_paintbrush_32x32Image = nav3_paintbrush_32x32.GetImage
getnav3_paintbrush_32x32Bitmap = nav3_paintbrush_32x32.GetBitmap
getnav3_paintbrush_32x32Icon = nav3_paintbrush_32x32.GetIcon

#----------------------------------------------------------------------
nav3_paintbrush_48x48 = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAADAAAAAwCAYAAABXAvmHAAAAAXNSR0IArs4c6QAAAARnQU1B"
    "AACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAAadEVYdFNvZnR3YXJlAFBhaW50Lk5F"
    "VCB2My41LjEwMPRyoQAAA4RJREFUaEPtmF1IU2EYx3e26TY33abOjzmFTczcdGGMjPwKbCRd"
    "SCjhIGleSEQOLFD6uBkyaAhSXdVFdWc3FvRxFUFB3URf0E1gdhfRB0TRFxiV9TviGZvndP8E"
    "+8PD83/Pu4u/r8/zvP9zTEUUUYQ8ZDIZU0tLS20gEBipqKiYrqurO19bW3uuvLx8enh42Lv+"
    "M3mIRqMVPp/vaFlZ2ROz2bzKoz/5YbPZTpHlIZ1OmxF+pLS09CPLAtFaOJ3OG2SFkAVO2uly"
    "uW5DDYWrwR/2PhKJVMFlYXJy0ky53IQaCteCHjhAlof6+vp9JEPRWpSUlDycmJgww+XB4XDc"
    "IxkK16KxsXEvWR7GxsYcVqt1BWooXA1q/0U8HrfC5aG/vz+kKIqhcC28Xu8UWRYQvZnJs218"
    "fNzFf+AHjwzF2+32pb6+PgdcDigJP5fRa+gVdc1sv0rSiac3ljs7O0NwOeCiqkTYc6gqUq39"
    "aCwWq2eMPlh/pt60n9xu93xvb6+btRx0dHQ4KYmc0PV4R3TOzs4qzPnWmpqajlQqVcozcbBw"
    "8qoNyBevxTfiNCbNT5YHRqVCSVyAGolfCxr6UWVlpTyHubi4aPJ4PFmooXA1KKv7TBqZ9riq"
    "qipF0lliLWje6+3t7WVweeACUj3OT8JQPP7mYjKZLIHLA29OOy0Wyz8vKHridFdXl0xzxqzf"
    "wu36GaoTrtoGyubY3NwcS4GgbIKIfwPViWfS/KKhD8Florm52UddL0GNxK/Q0CNwmWhoaHDi"
    "cTbesmtB2XyhJ3bDZWJwcNBGU96B6sRTTh+CwWAPXCba2toU6voyVCeesnnV1NQUgcsFdT1P"
    "0onn5JfC4bAsK7wReJfjJJ147PDjnp4emcZMg9/v30+J/IYWiGcK3e3u7vbA5SIUCu1hsuhe"
    "xF0u17WhoSEnXC5aW1u3c/K6Wxavf2l0dFSmr9EQiUQ24W/eQgvEc/uenZmZkelrNGB5/ZTN"
    "MjQnnPVqdXX1yWw2K+9jaz7i8bibEnkGzYlf9zWHFxYW1J/IRSKRsDHTb0HzT/47Uyih7ovG"
    "1NSUhRftgm82jMmvzP9d6r5opNNpBfFnoDnxqq8JBAI71H3x4JRPUCo58eqXNBxneG1TOvD1"
    "QZo09y5rt9tfcvLBtc3/AZz0VpIm/mksFpPtazYik8kozPckY/LgwMCAzE8eRRRRRB5Mpr9R"
    "zkZYLFqrQwAAAABJRU5ErkJggg==")
icon_index.append('nav3_paintbrush_48x48')
icon_catalog['nav3_paintbrush_48x48'] = nav3_paintbrush_48x48
getnav3_paintbrush_48x48Data = nav3_paintbrush_48x48.GetData
getnav3_paintbrush_48x48Image = nav3_paintbrush_48x48.GetImage
getnav3_paintbrush_48x48Bitmap = nav3_paintbrush_48x48.GetBitmap
getnav3_paintbrush_48x48Icon = nav3_paintbrush_48x48.GetIcon
        
        


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
        
        
class DemoMaskPanel(MaskPanelToolbar2):
    """
    Displays a single image panel with the ability to draw a mask overlaying
    the image.  User can also pan and zoom in on the image and turn on cross
    hairs to track the mouse cursor during paint and general mouse movements 
    modes.
    
    """
    # Activate event messages
    _EVENT_DEBUG = True
    
    def __init__( self, parent, tab, statusbar, **kwargs ):
        # statusbar has to be here for NavigationToolbar3Wx to discover on init()
        self.statusbar = statusbar
        # initiate plotter
        sizer = MaskPanelToolbar2.__init__(  self, 
                                             parent, 
                                             vertOn=True, 
                                             horizOn=True, 
                                             lcolor='gold',
                                             lw=0.5,
                                             **kwargs )  
        self.tab       = tab
        self.top       = wx.GetApp().GetTopWindow()
        self.parent    = parent
        self.statusbar = statusbar

    def on_motion(self, xloc, yloc, iplot):
        value = self.data[iplot][0]['data'][int(xloc)][int(yloc)]
        self.top.statusbar.SetStatusText( " Value = %s" % (str(value), ), 0)
        self.top.statusbar.SetStatusText( " X,Y = %i,%i" % (round(xloc),round(yloc)) , 1)
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

    def on_paint_motion(self, xloc, yloc, iplot):
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
         
        self.nb = wx.Notebook(self, -1, style=wx.BK_BOTTOM)
         
        panel1 = wx.Panel(self.nb, -1)
         
        self.view = DemoMaskPanel(panel1, self, self.statusbar, naxes=1, data=data)

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
                     ("Set Paint Radius=1", "", self.on_radius1),
                     ("Set Paint Radius=2", "", self.on_radius2),
                     ("Set Paint Radius=5", "", self.on_radius5),
                     ("Set Paint Radius=10","", self.on_radius10),
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
        
        data = [[data1], ]
        #data = [[data1], [data2]]

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
        
        data = [[data1], ]
        #data = [[data1], [data2]]

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
        
        data = [[data1], ]
        #data = [[data1], [data2]]

        self.view.set_data(data, keep_norm=keep_norm)
        self.view.update(no_draw=True, keep_norm=keep_norm)
        self.view.canvas.draw()        


    def on_radius1(self, event):
        self.view.paint_radius = 1        

    def on_radius2(self, event):
        self.view.paint_radius = 2        

    def on_radius5(self, event):
        self.view.paint_radius = 5        

    def on_radius10(self, event):
        self.view.paint_radius = 10        

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
            if i != 0: a[m1-i,:] = y    # Symmetriscal
        return a




if __name__ == '__main__':

    app   = wx.App( False )
    frame = MyFrame( title='WxPython and Matplotlib', size=(600,600) )
    frame.Show()
    app.MainLoop()

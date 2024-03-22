"""
Expansion of matplotlib embed in wx example by John Bender and Edward 
Abraham, see http://www.scipy.org/Matplotlib_figure_in_a_wx_panel

This version, plot_panel_spectrum.py, is a derivative of plot_panel.py that
has the specific purpose of displaying a plot of 1D spectral data in a 
variety of ways including: one spectrum, multiple lines in one spectrum, etc.  


Brian J. Soher, Duke University, March, 2014
"""

# Python modules

import math

# 3rd party modules
import matplotlib
import matplotlib.cm as cm
import wx
import numpy as np

# If we set the backend unconditionally, we sometimes get an undesirable
# message. 

if matplotlib.get_backend() != "WXAgg":
    matplotlib.use('WXAgg')

from matplotlib.transforms import blended_transform_factory
from matplotlib.patches    import Rectangle
from matplotlib.lines      import Line2D

# Our modules



class ImagePanel(wx.Panel):
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
                               color=None, 
                               bgcolor="#ffffff",
                               dpi=None, 
                               zoom='none', 
                               widlev=False,
                               middle=False,
                               unlink=False,
                               do_zoom_select_event=False,
                               do_zoom_motion_event=False,
                               do_widlev_select_event=False,
                               do_widlev_motion_event=False,
                               do_middle_select_event=False,
                               do_middle_motion_event=False,
                               do_scroll_event=False,
                               xscale_bump=0.0,
                               yscale_bump=0.0,
                               widlev_rate=3.0,     # multiplier for mouse motion
                               props_zoom=None,
                               data=None,
                               colormap=cm.gray,        # or gist_gray to reverse
                               ceiling=None,
                               floor=None,
                               **kwargs):

        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure

        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        self.parent = parent
        self.unlink = unlink
        self.xscale_bump = xscale_bump
        self.yscale_bump = yscale_bump
        self.widlev_rate = widlev_rate
        
        self.ceiling    = [ceiling for i in range(naxes)]              
        self.floor      = [floor   for i in range(naxes)]              
        self.width      = [256.0   for i in range(naxes)]
        self.level      = [128.0   for i in range(naxes)]
        self.width_base = [256.0   for i in range(naxes)]
        self.level_base = [128.0   for i in range(naxes)]
        self.norm       = matplotlib.colors.Normalize(vmin=0.0, vmax=255.0)
        self.cmap       = [colormap for i in range(naxes)]
        
        # AxesImage id returned from imshow() command
        self.imageid = [None for i in range(naxes)] 
        
        # store image data with ceil/floor applied, but not normalized 0-255
        # - used to display values from a given voxel not in grayscale range        
        self.img_hard_limits  = [None for i in range(naxes)]  
        
        # store image data with ceil/floor applied AND normalize to 0-255 
        # using max() and min() image value. 
        # - Start point for applying interactive width/level via widlev plugin
        self.img_norm = [None for i in range(naxes)]  
        
        # Under GTK we need to track self's size to avoid a continuous flow
        # of size events. 
        self._platform_is_gtk = ("__WXGTK__" in wx.PlatformInfo)
        self._current_size = (-1, -1)

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi, frameon=True )
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        
        
        # here we create the required naxes, add them to the figure, but we
        # also keep a permanent reference to each axes so they can be added
        # or removed from the figure as the user requests 1-N axes be displayed
        self.axes   = []
        for i in range(naxes):
            self.axes.append(self.figure.add_subplot(naxes,1,i+1))
        self.naxes = naxes   
        self.all_axes = list(self.axes)

        # internal data setup 
        if not data or len(data) != naxes:
            data = self._default_data()
        self.set_data(data)

        self.figure.subplots_adjust(left=0.01,right=0.99,
                                    bottom=0.01,top=0.99,
                                    wspace=0.0,hspace=0.0)

        # for images we don't show x or y axis values
        for axis in self.all_axes:
            axis.set_facecolor(bgcolor)
            axis.xaxis.set_visible(False)
            axis.yaxis.set_visible(False)
            
        self.zoom   = []
        self.widlev = []
        self.middle = []
        
        self.set_color( color )
        self._set_size()
        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._on_idle)
        self.Bind(wx.EVT_SIZE, self._on_size)

        # ensure that properties for zoom and reference regions exist
        if not props_zoom:
            props_zoom = dict(alpha=0.2, facecolor='yellow')

        #----------------------------------------------------------------------
        # enable Zoom, Reference, Middle and Scroll functionality as required

        if zoom == 'box':
            if not unlink:
                self.zoom = ZoomBox(  self, self.axes,
                                      drawtype='box',
                                      useblit=True,
                                      button=1,
                                      do_zoom_select_event=do_zoom_select_event,
                                      do_zoom_motion_event=do_zoom_motion_event,
                                      spancoords='data',
                                      rectprops=props_zoom)
            else:
                for axes in self.axes:
                    self.zoom.append(ZoomBox(  self, [axes],
                                          drawtype='box',
                                          useblit=True,
                                          button=1,
                                          do_zoom_select_event=do_zoom_select_event,
                                          do_zoom_motion_event=do_zoom_motion_event,
                                          spancoords='data',
                                          rectprops=props_zoom))
        if widlev:
            if not unlink:
                self.widlev = WidLevEvents(self, self.axes, button=3,
                                      do_widlev_select_event=do_widlev_select_event,
                                      do_widlev_motion_event=do_widlev_motion_event)
            else:
                for axes in self.axes:
                    self.widlev.append( WidLevEvents(self, [axes], button=3,
                                          do_widlev_select_event=do_widlev_select_event,
                                          do_widlev_motion_event=do_widlev_motion_event))
        
        if middle:
            if not unlink:
                self.middle = MiddleEvents(self, self.axes,
                                      do_middle_select_event=do_middle_select_event,
                                      do_middle_motion_event=do_middle_motion_event)
            else:
                for axes in self.axes:
                    self.middle.append( MiddleEvents(self, [axes],
                                          do_middle_select_event=do_middle_select_event,
                                          do_middle_motion_event=do_middle_motion_event))
          
        self.do_motion_event = True  
        self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)
        
        self.do_scroll_event = do_scroll_event
        if self.do_scroll_event:
            self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)

        # initialize plots with initial data and format axes
        self.set_data(self.data)
        self.update()

       


    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================

    def _dprint(self, a_string):
        if self._EVENT_DEBUG:
            print(a_string)


    def _on_size( self, event ):
        if self._platform_is_gtk:
            # This is a workaround 
            current_x, current_y = self._current_size
            new_x, new_y = tuple(event.GetSize())

            if (abs(current_x - new_x) > 1) or (abs(current_y - new_y) > 1):
                self._resizeflag = True
            else:
                # Size has only changed by one pixel or less. I ignore it.
                event.Skip(False)
        else:
            self._resizeflag = True


    def _on_idle( self, evt ):
        if self._resizeflag:
            self._resizeflag = False
            self._set_size()


    def _set_size( self ):
        pixels = tuple( self.parent.GetClientSize() )
        self.SetSize( pixels )
        self.canvas.SetSize( pixels )
        self.figure.set_size_inches( float( pixels[0] )/self.figure.get_dpi(),
                                     float( pixels[1] )/self.figure.get_dpi() )
        self._current_size = pixels


    def _on_move(self, event):
        """
        This is the internal method that organizes the data that is sent to the
        external user defined event handler for motion events. In here we 
        gather data values from line plots, determine 
        which axis we are in, then call the (hopefully) overloaded on_motion()
        method
        
        """
        if event.inaxes == None or not self.do_motion_event: return
        x0, y0, x1, y1 = bounds = event.inaxes.dataLim.bounds
        
        values, raw = self._get_values(event)
        iaxis = self._get_current_axis_index(event)
        self.on_motion(event.xdata, event.ydata, raw, bounds, iaxis)
        
        
    def _on_scroll(self, event):
        """
        This is the internal method that organizes the data that is sent to the
        external user defined event handler for scroll events. In here we 
        determine which axis we are in, then call the (hopefully) overloaded 
        on_scroll() method
        
        """
        if event.inaxes == None: return
        iaxis = self._get_current_axis_index(event)
        self.on_scroll(event.button, event.step, iaxis)        


    def _get_current_axis_index(self, event):
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i
        return iaxis


    def _default_data(self):
        data = []
        for i in range(self.naxes):
            data.append([self._dist(128),])
        return data


    def _get_values(self, event):
        """
        Generic utility function that polls the axes that the mouse is within
        to return a list of data values at the x location of the cursor.
        
        Note. value list here contains the uint8 scaled values in the 
          axes.images array while the rvalue list contains values from the
          data array after floor/ceil values are applied, but before the 
          window/level parameters are applied. Since there may be one or more
          data arrays in the canvas (hopefully with alpha values < 1.0) we
          return value and rvalue as a list. 
        
        """
        iaxis = self._get_current_axis_index(event)
        if iaxis is None: iaxis = 0
        
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
        xindx, yindx = self._get_data_index(event)

        value = []
        for image in event.inaxes.get_images():
            dat = image.get_array()
            if xindx < dat.shape[0] and yindx < dat.shape[1]:
                value.append(dat[xindx,yindx])
                    
        rvalue = []
        image = self.img_hard_limits[iaxis]
        if xindx < image.shape[0] and yindx < image.shape[1]:
            rvalue.append(image[xindx,yindx])
            
        return value, rvalue


    def _get_data_index(self, event, xvalue=None, yvalue=None):
        """
        Generic utility function that polls the axes that the mouse is within
        to return an index within the data array for the x location of the cursor.
        
        """
        indx = 0
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
        xpts,ypts = event.inaxes.get_images()[0].get_size()
    
        if not xvalue:
            xvalue = event.xdata
        if not yvalue:
            yvalue = event.ydata
    
        if xpts>=0:
            xindx = int(round((xpts-1) * (xvalue-x0)/x1))
        if ypts>=0:
            yindx = int(round((ypts-1) * (yvalue-y0)/y1))
    
        if xindx > (xpts-1): xindx = xpts-1
        if xindx < 0:        xindx = 0
    
        if yindx > (ypts-1): yindx = ypts-1
        if yindx < 0:        yindx = 0
    
        return xindx, yindx


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


    #=======================================================
    #
    #           User Accessible Data Functions  
    #
    #=======================================================

    def set_data(self, data, index=None):
        """
        User can set data into one or all axes using this method.
        
        If index is supplied, we assume that only one ndarray is being
        passed in via the data parameter. If no index is supplied then we
        assume that a list of ndarrays the size of self.naxes is being 
        passed in to replace all data in all axes.

        Example 1 - Data is a list of dicts
        
            raw  = {'data'  : raw_data,        # 2D numpy array 
                    'alpha' : 1.0 }            # value 0.0-1.0

            fit  = {'data' : fit_data,        # 2D numpy array 
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
                    if 'alpha' not in list(dat.keys()):
                        dat['alpha'] = 1.0
                    if 'cmap' not in list(dat.keys()):
                        dat['cmap'] = self.cmap[i]
                else:
                    # Only data in this item, so add all default values 
                    dat = { 'data'  : dat,
                            'alpha' : 1.0,
                            'cmap'  : self.cmap[i] }
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


    def update(self, index=None, keep_norm=False, no_draw=False):
        """
        Convenience function that runs through all the typical steps needed
        to refresh the screen after a set_data().
        
        The set_scale option is typically used only once to start set the 
        bounding box to reasonable bounds for when a zoom box zooms back 
        out.
        
        """
        self.apply_hard_limits(index=index)
        self.apply_norm_widlev(index=index, keep_norm=keep_norm)
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
            img_norm = self.img_norm[i]
            alpha    = ddict['alpha']
            cmap     = ddict['cmap']
                
            img_norm = self.calc_lut_value(img_norm, self.width[i], self.level[i])
                
            xmin, xwid, ymin, ywid = 0, img_norm.shape[1], 0, img_norm.shape[0]
            
            self.imageid[i] = axes.imshow(img_norm, norm=self.norm, 
                                                    cmap=cmap, 
                                                    alpha=alpha) 
            
            if xold != xwid or yold!=ywid or force_bounds:    
                # Set new bounds for dataLims to x,y extent of data in the 
                # new image. On reset zoom this is how far we reset the limits.
                axes.ignore_existing_data_limits = True
                axes.update_datalim([[xmin,ymin],[xmin+xwid,ymin+ywid]])
                # Matches viewLims view limits to the new data. By 
                # default, new data and updating causes display to show the 
                # entire image. Any zoom setting is lost.
                axes.set_xlim((xmin, xmin+xwid), auto=None)
                axes.set_ylim((ymin, ymin+ywid), auto=None)

            # may need this line in future if we do overlay
            #    self.figure.hold(True)
            #self.figure.hold(False)
            

    def calc_lut_value(self, data, width, level):
        """Apply Look-Up Table values to data for given width/level values."""

        conditions = [data <= (level-0.5-(width-1)/2), data > (level-0.5+(width-1)/2)]
        functions  = [0, 255, lambda data: ((data - (level - 0.5))/(width-1) + 0.5)*(255-0)]   # 3rd function is default
        lutvalue = np.piecewise(data, conditions, functions)

        # Convert the resultant array to an unsigned 8-bit array to create
        # an 8-bit grayscale LUT since the range is only from 0 to 255
        return np.array(lutvalue, dtype=np.uint8)


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


    def set_floor_ceiling(self, index, floor, ceiling):
        """
        Intent is to allow users to set floor and ceiling values for each image
        This can be abused in a variety of ways. Typically I would expect them 
        to use one of the scenarios below:
        
        1) provide one scalar index and a scalar floor and/or scalar ceiling
           value to be inserted
        2) set index to None and have floor and ceiling lists of naxes length 
           which map one to one to the existing arrays
        3) set index to None and have scalar floor and ceiling values that 
           map to all existing arrays.
        
        """
        indices = self.parse_indices(index)
        
        for i in indices:

            if ceiling is not None:
                if isinstance(ceiling, list):
                    if len(ceiling) == len(self.ceiling):
                        # apply one to one map of values to axes
                        self.ceiling[i] = ceiling[i]
                    else:
                        raise ValueError("Ceiling list wrong size. Set either all values or only one.")
                else:
                    # apply scalar value to all axes
                    self.ceiling[i] = ceiling 
    
            if floor is not None:
                if isinstance(floor, list):
                    if len(floor) == len(self.floor):
                        # apply one to one map of values to axes
                        self.floor[i] = floor[i]
                    else:
                        raise ValueError("Floor list wrong size. Set either all values or only one.")
                else:
                    # apply scalar value to all axes in indices
                    self.floor[i] = floor 
                            

    def apply_hard_limits(self, index=None):
        """ 
        Apply the user provided ceiling/floor values to limit the range of data
        values in the original image data.

        We save these floor to ceiling normalized images to serve as a starting 
        point from which we calculate normalized images. We also use the saved
        limited images to snag values out of to report in various events. 
               
        """
        indices = self.parse_indices(index)
        
        # loop through the data, apply ceiling/floor value and save new image
        for i in indices:

            data = self.data[i][0]['data'].copy()        
            orig_shape = data.shape 
            data = np.abs(data.flatten())
            if self.ceiling[i]: data[data>self.ceiling[i]] = self.ceiling[i]
            if self.floor[i]:   data[data<self.floor[i]]   = self.floor[i]
            data.shape = orig_shape
            self.img_hard_limits[i] = data       # ceil/floor applied
  
  
    def calc_norm_widlev(self, index=None):
        """ 
        Calculate the width and level values to that describe the original 
        data in images in terms of their data.min and data.max values

        """
        indices = self.parse_indices(index)
        for i in indices:
            data = self.img_hard_limits[i].copy()
            self.width_base[i] = int(np.abs(data.max()) + np.abs(data.min()))
            self.level_base[i] = int(self.width_base[i]*0.5 - np.abs(data.min()))


    def apply_norm_widlev(self, index=None, keep_norm=False):
        """ 
        Apply the calculated width and level values that describe the original 
        data in images in terms of their data.min and data.max values

        We save these 0-255 normalized images to serve as a starting point from 
        which we apply the user selected width and level values to window the 
        image on the screen.
               
        """
        indices = self.parse_indices(index)
        for i in indices:
            data = self.img_hard_limits[i].copy()
            if not keep_norm:
                self.calc_norm_widlev(i)
            self.img_norm[i] = self.calc_lut_value(data, self.width_base[i], self.level_base[i])


                       
    #=======================================================
    #
    #           User Accessible Plotting Functions  
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


    def change_naxes(self, n):
        """
        Allows user to determine serially which of the N axes are 
        included in the figure. Using this method the user supplies only the
        number of axes to include and the first 1:n axes in the long term
        storage list are added to the figure. 
        
        This method also updates the axes lists in the zoom, widlev and middle
        functor methods.
        
        """
        ncurrent = len(self.figure.axes)
        if n > self.naxes:
            return
        elif n < 0: 
            return
        elif n == ncurrent:
            return
        
        self.axes = self.all_axes[0:n]  

        # remove old axes, but don't destroy
        figure_axes = list(self.figure.axes)
        for axes in figure_axes:
            self.figure.delaxes(axes)
        
        # add back however many were requested
        for axes in self.axes:
            ax = self.figure.add_axes(axes)

        if not self.unlink:
            if self.zoom:
                self.zoom.axes = self.axes

            if self.widlev:
                self.widlev.axes = self.axes
                
            if self.middle:
                self.middle.axes = self.axes

        # this resets figure to have 1 or 2 or N axes shown
        naxes = len(self.axes)
        gs = matplotlib.gridspec.GridSpec(naxes,1)
        for i in range(naxes):
            self.figure.axes[i].set_position(gs[i].get_position(self.figure))
            self.figure.axes[i].set_subplotspec(gs[i])
            # bjs MPL deprecated self.figure.axes[i].change_geometry(naxes,1,i+1)


        self.canvas.draw()


    def display_naxes(self, flags):
        """
        Allows user to specify exactly which of the N axes defined in the 
        Init() method are included in the figure. 
        
        The user has to supply a boolean list of flags of the same length as
        the list of all_axes. The axes that correspond to flags set to True 
        are included in the figure.
        
        This method also updates the axes lists in the zoom, widlev and middle
        functor methods.
        
        """
        ncurrent = len(self.all_axes)
        nflags = len(flags)
        if nflags != ncurrent: return

        faxes = list(self.figure.axes)
        for axes in faxes:
            self.figure.delaxes(axes)
            
        for i, axes in enumerate(self.axes):
            if flags[i] != False:
                ax = self.figure.add_axes(axes)

        if not self.unlink:
            if self.zoom:
                self.zoom.axes = self.axes

            if self.widlev:
                self.widlev.axes = self.axes
                
            if self.middle:
                self.middle.axes = self.axes
        
        self.canvas.draw()


    def new_axes(self, axes):
        if isinstance(axes, list):
            self.axes = axes
        elif isinstance(axes, matplotlib.axes.Axes):
            self.axes = [axes]
        else:
            return

        if self.zoom is not None:
            self.zoom.new_axes(self.axes)

        if self.widlev is not None:
            self.widlev.new_axes(self.axes)

        if self.middle is not None:
            self.middle.new_axes(self.axes)
            
        if self.canvas is not self.axes[0].figure.canvas:
            self.canvas.mpl_disconnect(self.motion_id)
            self.canvas = self.axes[0].figure.canvas
            self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)       
        if self.figure is not self.axes[0].figure:
            self.figure = self.axes[0].figure


        
    #=======================================================
    #
    #           Default Event Handlers  
    #
    #=======================================================
        
    def on_motion(self, xdata, ydata, value, bounds, iaxis):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_move, xdata='+str(xdata)+'  ydata='+str(ydata)+'  val='+str(value)+'  bounds = '+str(bounds)+'  iaxis='+str(iaxis))
        
    def on_scroll(self, button, step, iaxis):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_move, button='+str(button)+'  step='+str(step)+'  iaxis='+str(iaxis))
        
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None, xdata=None, ydata=None ):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_select, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))
        
    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_move, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))

    def on_widlev_select(self, xstr, ystr, xend, yend, indx, reset=False):
        """ placeholder, overload for user defined event handling """
        self._dprint('ext on_widlev_select, X(str,end)='+str(xstr)+','+str(xend)+'  Y(str,end)='+str(ystr)+','+str(yend)+'  Index = '+str(indx))

    def on_widlev_motion(self, xcur, ycur, xprev, yprev, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_widlev_move, X(cur,prev)='+str(xcur)+','+str(xprev)+'  Y(cur,prev)='+str(ycur)+','+str(yprev)+'  Index = '+str(indx))

    def on_widlev_press(self, xloc, yloc, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_widlev_press, Xloc='+str(xloc)+'  Yloc='+str(yloc)+'  Index = '+str(indx))

    def on_middle_select(self, xstr, ystr, xend, yend, indx, reset=False):
        """ placeholder, overload for user defined event handling """
        self._dprint('ext on_middle_select, X(str,end)='+str(xstr)+','+str(xend)+'  Y(str,end)='+str(ystr)+','+str(yend)+'  Index = '+str(indx))
        
    def on_middle_motion(self, xcur, ycur, xprev, yprev, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_move, X(cur,prev)='+str(xcur)+','+str(xprev)+'  Y(cur,prev)='+str(ycur)+','+str(yprev)+'  Index = '+str(indx))

    def on_middle_press(self, xloc, yloc, indx, bounds=None, xdata=None, ydata=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_press, Xloc='+str(xloc)+'  Yloc='+str(yloc)+'  Index = '+str(indx))
        




class ZoomBox:
    """
    Select a min/max range of the x axes for a matplotlib Axes

    Example usage::

        from matplotlib.widgets import  RectangleSelector
        from pylab import *

        def onselect(xmin, xmax, value, ymin, ymax):
          'eclick and erelease are matplotlib events at press and release'
          print ' x,y min position : (%f, %f)' % (xmin, ymin)
          print ' x,y max position   : (%f, %f)' % (xmax, ymax)
          print ' used button   : ', eclick.button

        def toggle_selector(event):
            print ' Key pressed.'
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print ' RectangleSelector deactivated.'
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print ' RectangleSelector activated.'
                toggle_selector.RS.set_active(True)

        x = arange(100)/(99.0)
        y = sin(x)
        fig = figure
        axes = subplot(111)
        axes.plot(x,y)

        toggle_selector.RS = ZoomBox(axes, onselect, drawtype='line')
        connect('key_press_event', toggle_selector)
        show()
    """
    def __init__(self, parent, 
                       axes, 
                       drawtype='box',
                       minspanx=None, 
                       minspany=None, 
                       useblit=False,
                       lineprops=None, 
                       rectprops=None,
                       do_zoom_select_event=False, 
                       do_zoom_motion_event=False,
                       spancoords='data',
                       button=None):
        """
        Create a selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(pos_1, pos_2)

        and clear the drawn box/line. There pos_i are arrays of length 2
        containing the x- and y-coordinate.

        If minspanx is not None then events smaller than minspanx
        in x direction are ignored(it's the same for y).

        The rect is drawn with rectprops; default
          rectprops = dict(facecolor='red', edgecolor = 'black',
                           alpha=0.5, fill=False)

        The line is drawn with lineprops; default
          lineprops = dict(color='black', linestyle='-',
                           linewidth = 2, alpha=0.5)

        Use type if you want the mouse to draw a line, a box or nothing
        between click and actual position ny setting

        drawtype = 'line', drawtype='box' or drawtype = 'none'.

        spancoords is one of 'data' or 'pixels'.  If 'data', minspanx
        and minspanx will be interpreted in the same coordinates as
        the x and y axis, if 'pixels', they are in pixels

        button is a list of integers indicating which mouse buttons should
        be used for rectangle selection.  You can also specify a single
        integer if only a single button is desired.  Default is None, which
        does not limit which button can be used.
        Note, typically:
         1 = left mouse button
         2 = center mouse button (scroll wheel)
         3 = right mouse button
         
        """
        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []
        
        self.active = True                    # for activation / deactivation
        self.to_draw = []
        self.background = None
        self.axes_index = None

        self.do_zoom_select_event = do_zoom_select_event
        self.do_zoom_motion_event = do_zoom_motion_event
        
        self.useblit = useblit
        self.minspanx = minspanx
        self.minspany = minspany

        if button is None or isinstance(button, list):
            self.validButtons = button
        elif isinstance(button, int):
            self.validButtons = [button]

        assert(spancoords in ('data', 'pixels'))

        self.spancoords = spancoords
        self.eventpress = None          # will save the data (position at mouseclick)
        self.eventrelease = None        # will save the data (pos. at mouserelease)

        self.new_axes(axes, rectprops)


    def new_axes(self,axes, rectprops=None):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas
            
            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            
        if rectprops is None:
            rectprops = dict(facecolor='white', 
                             edgecolor= 'black',
                             alpha=0.5, 
                             fill=False)
        self.rectprops = rectprops

        for axes in self.axes:
            self.to_draw.append(Rectangle((0,0), 0, 1,visible=False,**self.rectprops))

        for axes,to_draw in zip(self.axes, self.to_draw):
            axes.add_patch(to_draw)
    

    def update_background(self, event):
        'force an update of the background'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        'return True if event should be ignored'
        # If ZoomBox is not active :
        if not self.active:
            return True

        # If canvas was locked
        if not self.canvas.widgetlock.available(self):
            return True

        # Only do selection if event was triggered with a desired button
        if self.validButtons is not None:
            if not event.button in self.validButtons:
                return True

        # If no button pressed yet or if it was out of the axes, ignore
        if self.eventpress == None:
            return event.inaxes not in self.axes

        # If a button pressed, check if the release-button is the same
        return  (event.inaxes not in self.axes or
                 event.button != self.eventpress.button)

    def press(self, event):
        'on button press event'
        # Is the correct button pressed within the correct axes?
        if self.ignore(event): return
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = False

        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i

        # make the drawn box/line visible get the click-coordinates,
        # button, ...
        for to_draw in self.to_draw:
            to_draw.set_visible(self.visible)
        self.eventpress = event
        return False


    def release(self, event):
        'on button release event'
        if self.eventpress is None or self.ignore(event): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = True
        
        # make the box/line invisible again
        for to_draw in self.to_draw:
            to_draw.set_visible(False)
        
        # left-click in place resets the x-axis or y-axis
        if self.eventpress.xdata == event.xdata and self.eventpress.ydata == event.ydata:
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds
            xdel = self.parent.xscale_bump*(x1-x0)
            ydel = self.parent.yscale_bump*(y1-y0)
            for axes in self.axes:
                axes.set_xlim(x0-xdel,x0+x1+xdel)
                axes.set_ylim(y0-ydel,y0+y1+ydel)
            self.canvas.draw()
            
            if self.do_zoom_select_event:
                self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True, 
                                                                                            iplot=self.axes_index, 
                                                                                            xdata=event.xdata,
                                                                                            ydata=event.ydata )
            
            return
        
        self.canvas.draw()
        # release coordinates, button, ...
        self.eventrelease = event

        if self.spancoords=='data':
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata

            # calculate dimensions of box or line get values in the right
            # order
        elif self.spancoords=='pixels':
            xmin, ymin = self.eventpress.x, self.eventpress.y
            xmax, ymax = self.eventrelease.x, self.eventrelease.y
        else:
            raise ValueError('spancoords must be "data" or "pixels"')

        # assure that min<max values
        if xmin>xmax: xmin, xmax = xmax, xmin
        if ymin>ymax: ymin, ymax = ymax, ymin
        # assure that x and y values are not equal
        if xmin == xmax: xmax = xmin*1.0001
        if ymin == ymax: ymax = ymin*1.0001

        spanx = xmax - xmin
        spany = ymax - ymin
        xproblems = self.minspanx is not None and spanx<self.minspanx
        yproblems = self.minspany is not None and spany<self.minspany
        if (xproblems or  yproblems):
            """Box too small"""    # check if drawed distance (if it exists) is
            return                 # not to small in neither x nor y-direction
        
        for axes in self.axes:
            axes.set_xlim((xmin,xmax))
            axes.set_ylim((ymin,ymax))
        self.canvas.draw()


        data_test = event.inaxes.images!=[]

        if self.do_zoom_select_event and data_test:
            # gather the values to report in a selection event
            value, raw = self.parent._get_values(event)            
            self.parent.on_zoom_select(xmin, xmax, raw, ymin, ymax, iplot=self.axes_index) # zeros are for consistency with box zoom

        self.axes_index   = None
        self.eventpress   = None              # reset the variables to their
        self.eventrelease = None              #   inital values
        
        return False


    def update(self):
        'draw using newfangled blit or oldfangled draw depending on useblit'
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, to_draw in zip(self.axes, self.to_draw):
                axes.draw_artist(to_draw)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()
        return False


    def onmove(self, event):
        'on motion notify event if box/line is wanted'
        if self.eventpress is None or self.ignore(event): return
        x,y = event.xdata, event.ydata              # actual position (with
                                                    #   (button still pressed)
        minx, maxx = self.eventpress.xdata, x       # click-x and actual mouse-x
        miny, maxy = self.eventpress.ydata, y       # click-y and actual mouse-y
        if minx>maxx: minx, maxx = maxx, minx       # get them in the right order
        if miny>maxy: miny, maxy = maxy, miny       
        for to_draw in self.to_draw:
            to_draw.set_x(minx)                    # set lower left of box
            to_draw.set_y(miny)
            to_draw.set_width(maxx-minx)           # set width and height of box
            to_draw.set_height(maxy-miny)

        data_test = event.inaxes.images!=[]
        
        if self.do_zoom_motion_event and data_test:
            # gather the values to report in a selection event
            value, raw = self.parent._get_values(event)
            self.parent.on_zoom_motion(minx, maxx, raw, miny, maxy, iplot=self.axes_index) # zeros are for consistency with box zoom
        
        self.update()
        return False

    def set_active(self, active):
        """ 
        Use this to activate / deactivate the RectangleSelector
        from your program with an boolean variable 'active'.
        
        """
        self.active = active

    def get_active(self):
        """ to get status of active mode (boolean variable)"""
        return self.active


class WidLevEvents:
    """
    Act on events having to do with image scaling width and level changes

    see below for example usage
    
    """

    def __init__(self, parent, axes, 
                       do_widlev_select_event=False, 
                       do_widlev_motion_event=False,
                       do_widlev_press_event=False,
                       button=None):
        """
        Attach this to a button and up/down motions while the button is
        pressed will change the level value in the ImagePanel,and left/right
        motions will change the width value in the ImagePanel

        button is a list of integers indicating which mouse buttons should
        be used for width and level changes.  You can also specify a single
        integer if only a single button is desired.  Default is None, which
        does not limit which button can be used.
        Note, typically:
         1 = left mouse button
         2 = center mouse button (scroll wheel)
         3 = right mouse button

        """
        self.parent = parent
        self.axes = None
        self.canvas = None
        self.cids = []

        self.background = None
        self.pressxy    = None
        self.axes_index = None

        self.do_widlev_select_event = do_widlev_select_event
        self.do_widlev_motion_event = do_widlev_motion_event
        self.do_widlev_press_event  = do_widlev_press_event

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prevxy = (0,0)

        if button is None or isinstance(button, list):
            self.validButtons = button
        elif isinstance(button, int):
            self.validButtons = [button]

        self.eventpress = None          # will save the data (position at mouseclick)

        self.new_axes(axes)


    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event',  self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event',   self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
    
    def ignore(self, event):
        'return True if event should be ignored'
        # Only do selection if event was triggered with a desired button
        if self.validButtons is not None:
            if not event.button in self.validButtons:
                return True

        # If no button pressed yet or if it was out of the axes, ignore
        if self.eventpress == None:
            return event.inaxes not in self.axes

        # If a button pressed, check if the release-button is the same
        return  (event.inaxes not in self.axes or
                 event.button != self.eventpress.button)

    def press(self, event):
        'on button press event'
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_widlev_motion_event:
            self.parent.do_widlev_event = False

        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i

        self.pressxy = event.x, event.y
        self.prevxy  = event.x, event.y

        if self.do_widlev_press_event:
            self.parent.on_widlev_press(event.x, event.y, self.axes_index)
        
        self.eventpress = event
        
        return False

    def release(self, event):
        'on button release event'
        if self.pressxy is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_widlev_motion_event:
            self.parent.do_widlev_event = True

        xstr, ystr = self.pressxy
        xend = event.x 
        yend = event.y 

        # left-click in place resets autoscale width/level values
        if self.eventpress.xdata == event.xdata and self.eventpress.ydata == event.ydata:
            indx = self.axes_index
            self.parent.width[indx] = 255.0
            self.parent.level[indx] = 128.0
            self.parent.update_images(index=indx)
            self.parent.canvas.draw()
            
            if self.do_widlev_select_event:
                self.parent.on_widlev_select(xstr, ystr, xend, yend, self.axes_index, True)
            return

        if self.do_widlev_select_event:
            self.parent.on_widlev_select(xstr, ystr, xend, yend, self.axes_index, False)

        self.axes_index = None
        self.pressxy    = None
        self.eventpress = None
        
        return False

    def onmove(self, event):
        'on motion notify event'
        if self.pressxy is None or self.ignore(event): return
        rate = self.parent.widlev_rate
        xcurr, ycurr = event.x, event.y
        xprev, yprev = self.prevxy
        
        self.prevxy = event.x, event.y

        indx  = self.axes_index
        xdelt = int((xprev-xcurr))  # user set divide to slow down change 
        ydelt = int((yprev-ycurr))  #  due to mouse motion.
        
        if abs(ydelt) >= abs(xdelt):    
            self.parent.level[indx] += ydelt
        else:                           
            self.parent.width[indx] += xdelt

        self.parent.update_images(index=indx)
#         for image in self.parent.img_hard_limits[indx]:
#             data_lut = self.parent.calc_lut_value(image, self.parent.width[indx], self.parent.level[indx])
#             self.parent.imageid[indx].set_data(data_lut)
#             self.parent.figure.hold(True)
#         self.parent.figure.hold(False)
        self.parent.canvas.draw()

        if self.do_widlev_motion_event:
            self.parent.on_widlev_motion(xcurr, ycurr, xprev, yprev, self.axes_index) 

        return False


class MiddleEvents:
    """
    Report events having to do with the middle button

    Example usage:

      axes = subplot(111)
      axes.plot(x,y)

      def onselect(vmin, vmax):
          print vmin, vmax
      middle = MiddleEvents(axes, onselect)

      onmove_callback is an optional callback that will be called on mouse move
      with the span range

    """
    def __init__(self, parent, axes, 
                               do_middle_select_event=False, 
                               do_middle_motion_event=False,
                               do_middle_press_event=False):
        """
        Create a span selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(vmin, vmax)

        and clear the span.

        If minspan is not None, ignore events smaller than minspan

        The span rect is drawn with rectprops; default
          rectprops = dict(facecolor='red', alpha=0.5)

        set the visible attribute to False if you want to turn off
        the functionality of the span selector

        """
        self.parent = parent
        self.axes = None
        self.canvas = None
        self.cids = []

        self.background = None
        self.pressxy = None
        self.axes_index = None

        self.do_middle_select_event = do_middle_select_event
        self.do_middle_motion_event = do_middle_motion_event
        self.do_middle_press_event  = do_middle_press_event

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prevxy = (0,0)

        self.new_axes(axes)


    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
    
    def ignore(self, event):
        'return True if event should be ignored'
        return  event.inaxes not in self.axes or event.button !=2

    def press(self, event):
        'on button press event'
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = False

        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i

        bounds = event.inaxes.dataLim.bounds

        self.pressxy = event.x, event.y
        self.prevxy  = event.x, event.y

        if self.do_middle_press_event:
            self.parent.on_middle_press(event.x, event.y, self.axes_index, bounds=bounds, xdata=event.xdata, ydata=event.ydata)
        
        self.eventpress = event
        
        return False

    def release(self, event):
        'on button release event'
        if self.pressxy is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = True

        xstr, ystr = self.pressxy
        xend = event.x 
        yend = event.y 

        # left-click in place resets the x-axis or y-axis
        if self.eventpress.xdata == event.xdata and self.eventpress.ydata == event.ydata:
            if self.do_middle_select_event:
                self.parent.on_middle_select(xstr, ystr, xend, yend, self.axes_index, True)
            
            return

        if self.do_middle_select_event:
            self.parent.on_middle_select(xstr, ystr, xend, yend, self.axes_index, False)

        self.axes_index = None
        self.pressxy = None
        self.eventpress = None
        
        return False

    def onmove(self, event):
        'on motion notify event'
        if self.pressxy is None or self.ignore(event): return
        xcurrent, ycurrent = event.x, event.y
        xprevious, yprevious = self.prevxy
        
        self.prevxy = event.x, event.y

        if self.do_middle_motion_event:
            self.parent.on_middle_motion(xcurrent, ycurrent, xprevious, yprevious, self.axes_index) 

        return False
    



#------------------------------------------------
# Test Code
#------------------------------------------------


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


class DemoImagePanel(ImagePanel):
    """Plots several lines in distinct colors."""

    # Activate event messages
    _EVENT_DEBUG = True
    
    def __init__( self, parent, tab, **kwargs ):
        # initiate plotter
        ImagePanel.__init__( self, parent, **kwargs )  
        self.tab    = tab
        self.top    = wx.GetApp().GetTopWindow()
        self.parent = parent
        self.count = 0

    def on_motion(self, xdata, ydata, value, bounds, iaxis):
        self.top.statusbar.SetStatusText( " Value = %s" % (str(value), ), 0)
        self.top.statusbar.SetStatusText( " X,Y = %i,%i" % (xdata,ydata) , 1)

    def on_scroll(self, button, step, iaxis):
        pass
        #self.set_vertical_scale(step)

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        delta  = xmax - xmin
        self.top.statusbar.SetStatusText(( " Point Range = %.2f to %.2f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText(( " dPoints = %i " % (delta, )), 2)

    def on_widlev_motion(self, xcur, ycur, xprev, yprev, indx):
        pass
        
    def on_widlev_select(self, xstr, ystr, xend, yend, indx, reset=False):
        pass
        



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
        self.view = DemoImagePanel( panel1, self, 
                                      naxes=2,
                                      zoom='box',
                                      widlev=True, 
                                      middle=False,
                                      unlink=False,
                                      do_zoom_select_event=True,
                                      do_zoom_motion_event=True,
                                      do_widlev_select_event=True,
                                      do_widlev_motion_event=True,
                                      do_middle_select_event=False,
                                      do_middle_motion_event=False,
                                      do_scroll_event=True,
                                      xscale_bump=0.0,
                                      yscale_bump=0.0,
                                      data = data,
                                      colormap=cm.gray,
                                      )
    
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        panel1.SetSizer(sizer)
        self.view.Fit()    
    
        self.nb.AddPage(panel1, "One")
    
        self.view.set_color( (255,255,255) )


    def menuData(self):
        return [("&File", (
                    ("", "", ""),
                    ("&Quit",    "Quit the program",  self.on_close))),
                 ("Tests", (
                     ("Show Both Images", "", self.on_show_two, wx.ITEM_RADIO),
                     ("Show only One",  "", self.on_show_one,   wx.ITEM_RADIO),
                     ("", "", ""),
                     ("Set Small Images - keep norm",  "", self.on_small_images_keep_norm),
                     ("Set Medium Images - keep norm", "", self.on_medium_images_keep_norm),
                     ("Set Large Images - keep norm",  "", self.on_large_images_keep_norm),
                     ("", "", ""),
                     ("Set Small Images",  "", self.on_small_images),
                     ("Set Medium Images", "", self.on_medium_images),
                     ("Set Large Images",  "", self.on_large_images),
                     ("", "", ""),
                     ("Recalc Norm",  "", self.on_recalc_norm),
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

    def on_show_one(self, event):
        self.view.change_naxes(1)

    def on_show_two(self, event):
        self.view.change_naxes(2)

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

        self.view.set_data(data)
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

        self.view.set_data(data)
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

        self.view.set_data(data)
        self.view.update(no_draw=True, keep_norm=keep_norm)
        self.view.canvas.draw()

    def on_recalc_norm(self, event):
        self.view.apply_norm_widlev(keep_norm=False)
        self.view.update_images()
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

#------------------------------------------------------------------------------

        
if __name__ == '__main__':

    app   = wx.App( 0 )
    frame = MyFrame( title='WxPython and Matplotlib', size=(600,600) )
    frame.Show()
    app.MainLoop()


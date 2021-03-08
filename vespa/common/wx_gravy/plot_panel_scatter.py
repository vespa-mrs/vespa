"""
Expansion of matplotlib embed in wx example by John Bender and Edward 
Abraham, see http://www.scipy.org/Matplotlib_figure_in_a_wx_panel

This version, plot_panel_scatter.py, is derived from PlotPanelPoints which
is a derivative of plot_panel.py and it has the purpose of displaying a plot 
of 1D x,y data in a variety of ways including: one line plot, a line plot 
with an overlay, multiple lines in one plot, etc. However, it was (is) not 
possible to make different markers in the same plot() different colors. The 
examples I found used the scatter() method to achieve this. Thus this version
of the plot_panel() package.

This version is more restrictive in that it requires a plot_option attibute
that describes the current state of the displayed data (e.g. whether x-axis
is shown, and if a zero line drawn, etc.).  These plot_option settings can be
controlled using the self.get_xxx() methods. This moves a lot of the 
functionality that had been handwritten over and over in various inheirited
classes into the base class.  

Of particular interest here is that we deal with real / imag / magn data types 
within the base class.  The user sets the plot_option.data_types attribute 
using the self.set_data_type_xxx methods and the base class does the rest. 

The data you want in your plot is passed in on __init__ and/or can be set anew
later using the self.set_data() method.  The data passed in is always a list 
of lists. However, each item in the list can vary as follows:
 list of list of ndarrays OR 
 list of list of dicts  OR
 list of list or ndarrays and dicts
 
 In any dict passed in, there has to be a 'data' entry that contains an ndarray
 You can also pass in the line colors to be applied for real, imag, and magn
 line display by setting entries for 'line_color_real', 'line_color_imaginary' and
 'line_color_magnitude'. If any of these are not present in a dict passed in, then 
 they are set using default values from the plot_options attribute.
 
 Each list of ndarrays (or dicts) can contain one or more arrays or dicts. When
 the data is displayed as Summed, each array is summed individually and then
 plotted. So if two ndarrays are passed in as part of one list, there will 
 always be at least two lines shown in each plot.  When data is not 'summed' 
 then all lines from all arrays are displayed.
 
 The only requirement for consistency across all ndarrays or dicts that are 
 passed in is that they all have the same number of plot points. Where
 the plot point dimension is that in the last entry of the ndarray.shape 
 value (ie. ndarray.shape[-1])

This version allows the user to zoom in on the figure using either a span selector 
or a box selector. You can also set a persistent span selector that acts as cursor 
references on top of whatever is plotted

ZoomSpan based on matplotlib.widgets.SpanSelector
CursorSpan based on matplotlib.widgets.SpanSelector
BoxZoom based on matplotlib.widgets.RectangleSelector

Brian J. Soher, Duke University, September, 2012
"""

# Python modules

import math
import warnings

# 3rd party modules
import matplotlib
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
from vespa.common.constants import DEGREES_TO_RADIANS, RADIANS_TO_DEGREES



def get_data_index(event, xvalue=None, reversex=False):
    """
    Generic utility function that polls the axes that the mouse is within
    to return an index within the data array for the x location of the cursor.
    
    """
    indx = 0
    x0, y0, x1, y1 = event.inaxes.dataLim.bounds
    npts = len(event.inaxes.lines[0].get_ydata())

    if not xvalue:
        xvalue = event.xdata

    if npts>=0 and (x1-x0)!=0:
        indx = int(round((npts-1) * (xvalue-x0)/(x1-x0)))
    else:
        indx = 0

    if reversex:        indx = npts - indx - 1
    if indx > (npts-1): indx = npts-1
    if indx < 0:        indx = 0

    return indx



class PlotPanelScatter(wx.Panel):
    """
    Derived from PlotPanelPoints. It was (is) not possible to make different
    markers in the same plot() different colors. The examples I found used
    the scatter() method to achieve this. 
    
    The PlotPanel has a Figure and a Canvas and 'n' Axes. The user defines
    the number of axes on Init and this number cannot be changed thereafter.
    However, the user can change the number of axes displaye in the Figure.
    
    Axes are specified on Init because the zoom and reference cursors 
    need an axes to attach to to init properly.  
    
    on_size events simply set a flag, and the actual resizing of the figure is 
    triggered by an Idle event.
    
    PlotPanel Functionality
    --------------------------------------------------
    left mouse - If zoom mode is 'span', click and drag zooms the figure.  
                 A span is selected along the x-axis. On release, the 
                 axes xlim is adjusted accordingly. If zoom mode is 'box', then
                 a zoom box is drawn during click and drag and figure is 
                 zoomed on both x- and y-axes upon release.
               
    left mouse - click in place, un-zooms the figure to maximum x-data or
                 x-data and y-data bounds.
                 
    right mouse - If reference mode is True/On, then click and drag will draw
                  a span selector in the canvas that persists after release.
                  
    middle mouse - (or scroll button click), if do_middle_select_event and/or 
                   do_middle_motion_event are True then select, release and 
                   motion events are returned for these uses of the middle
                   mouse button. Mouse location and axes index values are 
                   returned
                   
    scroll roll - if do_scroll_event is True then these events are returned if 
                  they occur within an axes. Mouse location and axes index
                  values are returned

    """

    # Set _EVENT_DEBUG to True to activate printing of messages to stdout 
    # during events.
    _EVENT_DEBUG = False

    def __init__(self, parent, naxes=2,
                               color=None, 
                               dpi=None, 
                               reversex=False,
                               zoom='none', 
                               reference=False, 
                               middle=False,
                               unlink=False,
                               zoom_button=1,
                               middle_button=2,
                               refs_button=3,
                               do_zoom_select_event=False,
                               do_zoom_motion_event=False,
                               do_refs_select_event=False,
                               do_refs_motion_event=False,
                               do_middle_select_event=False,
                               do_middle_motion_event=False,
                               do_middle_press_event=False,
                               do_scroll_event=False,
                               xscale_bump=0.0,
                               yscale_bump=0.0,
                               props_zoom=None,
                               props_cursor=None,
                               data=None,
                               prefs=None,
                               line_width=1.0,
                               **kwargs):

        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure

        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        if zoom:
            if reference:
                if refs_button==zoom_button:
                    raise ValueError('Zoom and Reference button numbers the same.')
            if middle:
                if middle_button==zoom_button:
                    raise ValueError('Zoom and Middle button numbers the same.')
        if reference:
            if middle:
                if middle_button==refs_button:
                    raise ValueError('Reference and Middle button numbers the same.')
        self.zoom_button = zoom_button
        self.refs_button = refs_button
        self.middle_button = middle_button

        self.parent = parent
        self.reversex = reversex
        self.unlink = unlink
        self.xscale_bump = xscale_bump
        self.yscale_bump = yscale_bump
        
        # Under GTK we need to track self's size to avoid a continuous flow
        # of size events. 
        self._platform_is_gtk = ("__WXGTK__" in wx.PlatformInfo)
        self._current_size = (-1, -1)

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi )
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        self.figure.set_tight_layout({'pad': 0.1, 'h_pad': 0.1, 'w_pad': 0.1})
        
        # here we create the required naxes, add them to the figure, but we
        # also keep a permanent reference to each axes so they can be added
        # or removed from the figure as the user requests 1-N axes be displayed
        self.axes   = []
        for i in range(naxes):
            self.axes.append(self.figure.add_subplot(naxes,1,i+1))
        self.naxes = naxes   
        self.all_axes = list(self.axes)

        # plot format setup
        if not prefs:
            prefs = fake_prefs()
        self.prefs = prefs

        # ensure that all prefs are present
        if not hasattr(self.prefs,"data_type_summed"):
            self.prefs.data_type_summed = False

        # internal data setup 
        if not data or len(data) != naxes:
            data = self._default_data()
        self.set_data(data)
        self.multiplier       = [1     for i in range(naxes)]
        self.ref_locations    = 0,1

        self.xlabel           = 'Points'
        self.vertical_scale   = 1.0
        self.dataymax         = 1.0

        data_type = 'real'
        if hasattr(prefs, 'data_type_magnitude'):
            if prefs.data_type_magnitude: data_type = 'magnitude'
        if hasattr(prefs, 'data_type_imaginary'):
            if prefs.data_type_imaginary: data_type = 'imaginary'
        if hasattr(prefs, 'data_type_real'):
            if prefs.data_type_real: data_type = 'real'

        self.data_type        = [data_type for i in range(naxes)]
        self.data_type_summed = [self.prefs.data_type_summed for i in range(naxes)]
        self.line_width       = [self.prefs.line_width for i in range(naxes)]

        # control display of weight array in plots
        self.overlay1            = np.zeros(self.dim0)
        self.show_overlay1       = [False for i in range(naxes)]

        for axis in self.all_axes:
            axis.set_facecolor(self.prefs.bgcolor)
            
        self.zoom = []
        self.refs = []
        self.middle = []
        
        self.set_color( color )
        self._set_size()
        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._on_idle)
        self.Bind(wx.EVT_SIZE, self._on_size)

        # ensure that properties for zoom and reference regions exist
        if not props_zoom:
            props_zoom = dict(alpha=0.2, facecolor='yellow')

        if not props_cursor:
            props_cursor = dict(alpha=0.2, facecolor='purple')

        #----------------------------------------------------------------------
        # enable Zoom, Reference, Middle and Scroll functionality as required

        if zoom == 'span':  
        
            if not unlink:
                self.zoom = ZoomSpan( self, self.all_axes,
                                      useblit=True,
                                      button=zoom_button,
                                      do_zoom_select_event=do_zoom_select_event,
                                      do_zoom_motion_event=do_zoom_motion_event,
                                      rectprops=props_zoom)
            else:
                for axes in self.axes:
                    self.zoom.append( ZoomSpan( self, [axes],
                                          useblit=True,
                                          button=zoom_button,
                                          do_zoom_select_event=do_zoom_select_event,
                                          do_zoom_motion_event=do_zoom_motion_event,
                                          rectprops=props_zoom))
        if zoom == 'box':
            if not unlink:
                self.zoom = ZoomBox(  self, self.axes,
                                      drawtype='box',
                                      useblit=True,
                                      button=zoom_button,
                                      do_zoom_select_event=do_zoom_select_event,
                                      do_zoom_motion_event=do_zoom_motion_event,
                                      spancoords='data',
                                      rectprops=props_zoom)
            else:
                for axes in self.axes:
                    self.zoom.append(ZoomBox(  self, [axes],
                                          drawtype='box',
                                          useblit=True,
                                          button=zoom_button,
                                          do_zoom_select_event=do_zoom_select_event,
                                          do_zoom_motion_event=do_zoom_motion_event,
                                          spancoords='data',
                                          rectprops=props_zoom))
        if reference:
            if not unlink:
                self.refs = CursorSpan(self, self.axes,
                                      useblit=True,
                                      button=refs_button,
                                      do_refs_select_event=do_refs_select_event,
                                      do_refs_motion_event=do_refs_motion_event,
                                      rectprops=props_cursor)
            else:
                for axes in self.axes:
                    self.refs.append(CursorSpan(self, [axes],
                                          useblit=True,
                                          button=refs_button,
                                          do_refs_select_event=do_refs_select_event,
                                          do_refs_motion_event=do_refs_motion_event,
                                          rectprops=props_cursor))
        if middle:
            if not unlink:
                self.middle = MiddleEvents(self, self.axes,
                                           button=middle_button,
                                           do_middle_press_event=do_middle_press_event,
                                           do_middle_select_event=do_middle_select_event,
                                           do_middle_motion_event=do_middle_motion_event)
            else:
                for axes in self.axes:
                    self.middle.append( MiddleEvents(self, [axes],
                                          button=middle_button,
                                          do_middle_select_event=do_middle_select_event,
                                          do_middle_motion_event=do_middle_motion_event,
                                          do_middle_press_event=do_middle_press_event))
        
        self.do_motion_event = True  
        self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)
        
        self.do_scroll_event = do_scroll_event
        if self.do_scroll_event:
            self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)

        self.shift_is_held = False      # flag for shift+ButtonLeft or shift+ButtonRight
        self.canvas.mpl_connect('key_press_event', self._on_key_press)
        self.canvas.mpl_connect('key_release_event', self._on_key_release)

        # initialize plots with initial data and format axes
        self.set_data(self.data)
        self.update(set_scale=True)
        

    @property
    def dim0(self):
        return self.data[0][0]['data'].shape[-1]

        


    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================
            
    def _calculate_scale(self):
        """ 
        This is usually a one time call to set up various ylim values on the
        plot_panel. Subsequently, the menu_events take care of setting
        these options and then do a canvas.plot() call to refresh        
        """

        # take min/max only from first data set, since it will always be there
        xx = self.all_axes[0].lines[0].get_xdata()
        xmin = min(xx)
        xmax = max(xx)

        if xmin == xmax:
            xmin, xmax = xmin-0.05, xmax+0.05
        ymax = max(np.abs(self.data[0][0]['data'].flatten()))
        ymin = -ymax

        if ymin == ymax == 0: ymax = 1.0
        self.dataymax = ymax
        self.vertical_scale = ymax

        for axes in self.all_axes:

            # ensure bounds are correct on start 
            axes.ignore_existing_data_limits = True
            if self.prefs.zero_line_plot_top:
                axes.update_datalim([[xmin,-self.dataymax],[xmax,self.dataymax*0.1]])
            elif self.prefs.zero_line_plot_middle:
                axes.update_datalim([[xmin,-self.dataymax],[xmax,self.dataymax]])
            elif self.prefs.zero_line_plot_bottom:
                axes.update_datalim([[xmin,-self.dataymax*0.1],[xmax,self.dataymax]])
    
            x0, y0, x1, y1 = axes.dataLim.bounds
            xdel = self.xscale_bump*(x1-x0)
            ydel = self.yscale_bump*(y1-y0)
            for axes in self.axes:
                if self.reversex:
                    axes.set_xlim(x0+x1+xdel,x0-xdel)
                else:
                    axes.set_xlim(x0-xdel,x0+x1+xdel)
                axes.set_ylim(y0-ydel,y0+y1+ydel)


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
        
        values = self.get_values(event, reversex=self.reversex)
        
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i

        self.on_motion(event.xdata, event.ydata, values, bounds, iaxis)
        
        
    def _on_scroll(self, event):
        """
        This is the internal method that organizes the data that is sent to the
        external user defined event handler for scroll events. In here we 
        determine which axis we are in, then call the (hopefully) overloaded 
        on_scroll() method
        
        """
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i
                
        self.on_scroll(event.button, event.step, iaxis)        


    def _on_key_press(self, event):
        """internal method to set up shift+button events"""
        if event.key == 'shift':
            self.shift_is_held = True

    def _on_key_release(self, event):
        """internal method to release shift+button events"""
        if event.key == 'shift':
            self.shift_is_held = False


    def _default_data(self):
        data = []
        for i in range(self.naxes):
            data.append([np.zeros([1,50],)])
        return data


    def get_values(self, event, reversex=False):
        """
        Generic utility function that polls the axes that the mouse is within
        to return a list of data values at the x location of the cursor.
        
        """
        value = []
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
        indx = get_data_index(event, reversex=reversex)
        data = event.inaxes.lines[0].get_ydata()
        if len(data.shape)==0:
            value = [data]
        else:
            for line in event.inaxes.lines:
                dat = line.get_ydata()
                if indx < len(dat):
                    value.append(dat[indx])
        
        if value == []: value = [0.0]
        return value


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
        
            raw  = {'data' : raw_data,                 # complex numpy array 
                    'line_color_real'      : 'blue',
                    'line_color_imaginary' : 'red',
                    'line_color_magnitude' : 'blue',
                    'markevery'            : [],
                    'markevery_color'      : 'purple',
                    'xaxis_values' : np.arange(len(raw_data))/1000.0 }    # set x-values in [msec]

            fit  = {'data' : fit_data, 
                    'line_color_real'      : 'black',
                    'line_color_imaginary' : 'purple',
                    'line_color_magnitude' : 'black',
                    'markevery'            : [],
                    'markevery_color'      : 'green'    }    # default xaxis values used here

            dif  = {'data' : fit_data - raw_data, 
                    'line_color_real'      : 'green',
                    'line_color_imaginary' : 'green',
                    'line_color_magnitude' : 'green' }    # default xaxis values used here

            data = [raw, fit, dif]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized, no_draw=True)
            self.view.canvas.draw()

        Example 2 - Data is a single numpy array, the colors dict will use
                    default values set in set_data() method
        
            data = [raw_data,]        # complex numpy array
            data = [[data]]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized, no_draw=True)
            self.view.canvas.draw()
                    
        ---------------------------------------------------------------------
                    
        from  https://stackoverflow.com/questions/32653825/how-to-color-data-points-based-on-some-rules-in-matplotlib?noredirect=1&lq=1
                    
        import matplotlib.pyplot as plt
        import numpy as np
        
        k=np.array([12,11,12,12,20,10,12,0,12,10,11])
        x2=np.arange(1,12,1)
        
        # find out which parameters are more than 1.5*std away from mean
        warning = np.abs(k-np.mean(k)) > 1.5*np.std(k)
        
        # enable drawing of multiple graphs on one plot
        plt.hold(True)
        
        # draw some lines behind the scatter plots (using zorder)
        plt.plot(x2, k, c='black', zorder=-1)
        
        # scatter valid (not warning) points in blue (c='b')
        plt.scatter(x2[~warning], k[~warning], label='signal', c='b')
        
        # scatter warning points in red (c='r')
        plt.scatter(x2[warning], k[warning], label='warning', c='r')
        
        # draw the legend
        plt.legend()
        
        # show the figure
        plt.show()                    
                    
                    
                    
        """
        for i, item in enumerate(data):
            for j, dat in enumerate(item):
                if isinstance(dat, dict):
                    # Dict in this item, but ensure all line color types have defaults
                    if 'line_color_real' not in list(dat.keys()):
                        dat['line_color_real'] = self.prefs.line_color_real
                    if 'line_color_imaginary' not in list(dat.keys()):
                        dat['line_color_imaginary'] = self.prefs.line_color_imaginary
                    if 'line_color_magnitude' not in list(dat.keys()):
                        dat['line_color_magnitude'] = self.prefs.line_color_magnitude
                    if 'markevery' not in list(dat.keys()):
                        dat['markevery'] = []
                    if 'markevery_color' not in list(dat.keys()):
                        dat['markevery_color'] = 'green'
                else:
                    # Only data in this item, so add all default line colors 
                    dat = {'data' : dat,
                           'line_color_real'      : self.prefs.line_color_real,
                           'line_color_imaginary' : self.prefs.line_color_imaginary,
                           'line_color_magnitude' : self.prefs.line_color_magnitude,
                           'markevery'            : [],
                           'markevery_color'      : 'green'}

                item[j] = dat
        
        
        if index:
            if index < 0 or index >= self.naxes:
                raise ValueError("index must be within that number of axes in the plot_panel")
            
            if len(data[0][0]['data'].shape[-1]) != self.data[0][0]['data'].shape[-1]:
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
                
                    if d.shape[-1] != data[0][0]['data'].shape[-1]:
                        raise ValueError("all ndarrays must have same number of spectral points")
            
            self.data = data

        # add default x-axis values if needed
        data = self.data
        for i, item in enumerate(data):
            for j, dat in enumerate(item):
                if 'xaxis_values' not in list(dat.keys()):
                    dat['xaxis_values'] = np.arange(len(dat['data'][0]))
                else:
                    # xaxis length must match data dimension
                    if len(dat['xaxis_values']) != len(dat['data'][0]):
                        dat['xaxis_values'] = np.arange(len(dat['data'][0]))
                item[j] = dat
        self.data = data



    def update(self, set_scale=False, no_draw=False):
        """
        Convenience function that runs through all the typical steps needed
        to refresh the screen after a set_data().
        
        The set_scale option is typically used only once to start set the 
        bounding box to reasonable bounds for when a zoom box zooms back 
        out.
        
        """
        self.update_plots()
        if set_scale:
            self._calculate_scale()
        self.update_axes()
        if not no_draw:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.canvas.draw()


    def update_axes(self):
        """
        Convenience function that runs through all the typical steps needed
        to refresh the screen after an axes change.  self.all_axes[0].lines[0].get_ydata()
        
        """
        self.format_axes()
        self.set_ylim()
    
    
    def update_plots(self):
        """
        Sets the data from the numpy arrays into the axes.
        
        Eventually, this will include a step to copy the data into a temp
        buffer where phase or other actions can be applied without messing
        up the original data.
        
        The only formatting performed here is to ensure that the x-axis
        has proper xmin and xmax set to follow the self.reversex setting
        
        """
        self.reversex = False

        #xx  = np.arange(self.dim0)
        
        # NB. we only use first set of x-axis values for now
        xx  = self.data[0][0]['xaxis_values']
        
        # take min/max only from first axes, since it will always be there
        ymax = max(np.abs(self.data[0][0]['data'].flatten()))
        ymin = -ymax
        xmax = max(xx)
        xmin = min(xx)

        if xmin == xmax:
            xmin, xmax = xmin-0.001, xmax+0.001

        for i, axes in enumerate(self.all_axes):

            # store current xlim values to restore later if in new range
            if self.reversex:
                old_xmax, old_xmin = axes.get_xlim()
            else:
                old_xmin, old_xmax = axes.get_xlim()

            axes.lines = []
            width = self.line_width[i]

            dlist = self.data[i]
            
            for ddict in dlist:
                
                data = ddict['data'].copy()
                if self.data_type_summed[i]:
                    data = np.sum(data, 0)
                    data.shape = 1,data.shape[0]
    
                for j, _ in enumerate(data):
                    color = ddict['line_color_real']
                    axes.plot(xx, data[j,:], color=color, linewidth=width)
                    if ddict['markevery']:
                        axes.plot(xx, data[j,:], linestyle='', marker='o', markevery=ddict['markevery'], color=ddict['markevery_color'])
                    
            # zero line 
            axes.axhline(0, color=self.prefs.zero_line_plot_color,
                            linestyle=self.prefs.zero_line_plot_style,
                            linewidth=width)  

            # overlay1 line 
            if len(self.overlay1) != data.shape[-1]:
                self.overlay1 = np.zeros(data.shape[-1])
            axes.plot(xx, self.overlay1, color=self.prefs.zero_line_plot_color,
                            linestyle=self.prefs.zero_line_plot_style,
                            linewidth=width)  

            # if x-axis has changed, ensure bounds are appropriate
            x0, y0, x1, y1 = axes.dataLim.bounds
            axes.ignore_existing_data_limits = True  
            axes.update_datalim([[xmin,y0],[xmax,y1+y0]])  
            if self.reversex:
                if old_xmin < xmin or old_xmax > xmax:
                    axes.set_xlim(xmax,xmin)
                else:
                    axes.set_xlim(old_xmax,old_xmin)
            else:
                if old_xmin < xmin or old_xmax > xmax:
                    axes.set_xlim(xmin,xmax)
                else:
                    axes.set_xlim(old_xmin,old_xmax)


    def get_data(self, index):
        """ Return a copy of one and only one of the data sets. """
        if index < 0 or index >= self.naxes:
            return None
        else:
            return self.data[index][0]['data'].copy()


    def format_axes(self):
        """
        Here the plot_option settings are enforced for display of x-axis
        label, data_type (real/imag/magn) selection, display of zero line
        and position of zero line.
        
        """
        xx = self.all_axes[0].lines[0].get_xdata()
        xmin = min(xx)
        xmax = max(xx)
        ymax = self.vertical_scale 
        ymin = -ymax

        if xmin == xmax:
            xmin, xmax = xmin-0.05, xmax+0.05
        
#        if self.prefs.xaxis_show:
#            self.figure.subplots_adjust(left=0.0,right=0.999,
#                                        bottom=0.15,top=1.0,
#                                        wspace=0.0,hspace=0.0)
#            self.all_axes[self.naxes-1].xaxis.set_visible(True)                                             
#        else:
#            self.figure.subplots_adjust(left=0.0,right=0.999,
#                                        bottom=0.0,top=1.0,
#                                        wspace=0.0,hspace=0.0)  
#            self.all_axes[self.naxes-1].xaxis.set_visible(False)
#
#        self.all_axes[self.naxes-1].set_xlabel(self.xlabel)

        # set axes on/off and use appropriate title
        naxes = len(self.axes) - 1
        for i,axes in enumerate(self.axes):
            xaxis_show = self.prefs.xaxis_show if i==naxes else False
            the_label = self.xlabel if i==naxes else ''
            axes.xaxis.set_visible(xaxis_show)
            axes.yaxis.set_visible(False)
            axes.set_xlabel(the_label)

        for j, axes in enumerate(self.all_axes):

            # flag whether to display overlay1 line
            axes.lines[int((len(axes.lines) - 1))].set_visible(self.show_overlay1[j])

            # flag whether to display zero line
            axes.lines[int((len(axes.lines) - 2))].set_visible(self.prefs.zero_line_plot_show)

            # set zero line at top/middle/bottom of plot
            axes.ignore_existing_data_limits = True
            if self.prefs.zero_line_plot_top:
                axes.update_datalim([[xmin,-self.dataymax],[xmax,self.dataymax*0.1]])
            elif self.prefs.zero_line_plot_middle:
                axes.update_datalim([[xmin,-self.dataymax],[xmax,self.dataymax]])
            elif self.prefs.zero_line_plot_bottom:
                axes.update_datalim([[xmin,-self.dataymax*0.1],[xmax,self.dataymax]])


    def set_ylim(self):
        """
        Refresh axes ylim values based on plot_option zero_line_location
        setting
        
        """
        if self.prefs.zero_line_plot_bottom:
            ymax =  1.0*self.vertical_scale
            ymin = -0.1*ymax
        elif self.prefs.zero_line_plot_middle: 
            ymax =  1.0*self.vertical_scale 
            ymin = -ymax 
        elif self.prefs.zero_line_plot_top: 
            ymin = -1.0*self.vertical_scale 
            ymax = -0.1*ymin 
        ylim = (ymin,ymax)
        if ymax == ymin:
            ylim = (1,0)
        for axes in self.all_axes:
            axes.set_ylim(ylim)


    def reset_xlim(self):
        """ set xlim values to max and min bounding box """
        for i, axes in enumerate(self.all_axes):  
            x0, y0, x1, y1 = axes.dataLim.bounds
            if self.reversex:
                axes.set_xlim(x0+x1,x0)
            else:
                axes.set_xlim(x0,x0+x1)

    def reset_ylim(self):
        """ set ylim values to max and min bounding box """
        for axes in self.all_axes:  
            x0, y0, x1, y1 = axes.dataLim.bounds
            axes.set_ylim(y0,y0+y1)

        
    def set_data_type_summed(self, index=None):
        if index is not None:
            for indx in index:
                if indx >=0 and indx < self.naxes:
                    self.data_type_summed[indx] = not self.data_type_summed[indx]
        else:
            self.data_type_summed = [not item for item in self.data_type_summed] 

        
    def set_vertical_scale(self, step, scale_mult=1.25):
        '''
        Vertical scale events are typically controlled by the scroll event
        that return +/- values for step depending on roll being up or
        down. Also, some computers are set to increment each step by x3 or 
        more (e.g. how many lines in Word or Browser to scroll with each 
        click).
        
        So we have to keep step from getting too large. To do this we set a 
        flag to indicate whether we are scrolling up or down and then set the 
        step value to be positive and divide by 3 but don't let any given step 
        be less than 1.
        
        '''
        shrink = False
        if step < 0: 
            shrink = True
            step = -step
        step = step/3
        if step < 1: step = 1
        
        if shrink:
            self.vertical_scale = self.vertical_scale * scale_mult * step
        else:
            self.vertical_scale = self.vertical_scale / (scale_mult * step)
              
        self.set_ylim()
        self.canvas.draw() 


    def set_overlay1(self, value):
        if len(value) == self.dim0:
            self.overlay1 = value


    def set_show_overlay1(self, value, index=-1):
        """Given a boolean and an optional index, sets the overlay1 flag to
        the boolean value for the given axis. If no index is given, the flag
        is set for all axes.
        """
        if index == -1:
            # Set all axes
            self.show_overlay1 = [value for i in self.all_axes]
        else:
            self.show_overlay1[index] = value


    def calculate_area(self):
        ''' 
        Calculates & returns the area and rms area selected between the 
        reference lines.

        '''
        rstr, rend = self.ref_locations 
        if rstr > rend: rstr, rend = rend, rstr

        rstr = int(rstr)
        rend = int(rend)

        all_areas = []
        all_rms = []
        for axes in self.all_axes:
            if axes.lines != []:
                data = axes.lines[0].get_ydata()
                all_areas.append(sum(data[rstr:rend])) 
                all_rms.append(self.rms(data[rstr:rend]))
            else:
                all_areas.append(0.0)
                all_rms.append(0.0)
            
        return all_areas, all_rms       
    
    
    def rms(self, data):
        if len(data) < 3:
            return 0.0
        return np.sqrt(np.sum( (data - np.mean(data))**2 )/(len(data)-1.0) )  
                        
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

    def refresh_cursors(self):
        """ redraws the reference cursor span on user request """
        if self.refs == None: return
        if self.refs.rect == []: return
        for axes, rect in zip(self.axes, self.refs.rect):
            axes.add_patch(rect)
        self.canvas.draw()

    def change_naxes(self, n):
        """
        Allow user to determine interactively which of the N axes are included
        in the figure. NB. This method irrevocably removes axes from the figure.
        User supplies only the number of axes to include and the first 1:n
        axes in the long term storage list are retained in the figure. This method
        also updates the axes lists in any zoom, refs or middle objects.
        
        """
        if n > self.naxes or n < 0 or n == len(self.figure.axes):
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
                
            if self.refs:
                self.refs.axes = self.axes
    
            if self.middle:
                self.middle.axes = self.axes

        # this resets figure to have 1 or 2 or N axes shown
        naxes = len(self.axes)
        for i in range(naxes):
            self.figure.axes[i].change_geometry(naxes,1,i+1)

        self.canvas.draw()


    def display_naxes(self, flags):
        """
        Allows user to specifiy exactly which of the N axes defined in the 
        init() method are included in the figure. 
        
        The user has to supply a boolean list of flags of the same length as
        the list of all_axes. The axes that correspond to flags set to True 
        are included in the figure.
        
        This method also updates the axes lists in the zoom, refs and middle
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
                
            if self.refs:
                self.refs.axes = self.axes
    
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
        if self.reference is not None:
            self.refs.new_axes(self.axes)
            
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
        
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_select, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))
        
    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_move, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))

    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None, rindex=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_refs_select, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val))
        
    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_refs_move, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val))

    def on_middle_select(self, xstr, ystr, xend, yend, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('ext on_middle_select, X(str,end)='+str(xstr)+','+str(xend)+'  Y(str,end)='+str(ystr)+','+str(yend)+'  Index = '+str(indx))
        
    def on_middle_motion(self, xcur, ycur, xprev, yprev, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_move, X(cur,prev)='+str(xcur)+','+str(xprev)+'  Y(cur,prev)='+str(ycur)+','+str(yprev)+'  Index = '+str(indx))

    def on_middle_press(self, xloc, yloc, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_press, Xloc='+str(xloc)+'  Yloc='+str(yloc)+'  Index = '+str(indx))
        





class ZoomSpan:
    """
    Select a min/max range of the x or y axes for a matplotlib Axes

    Example usage:

      axes = subplot(111)
      axes.plot(x,y)

      def onselect(vmin, vmax):
          print(vmin, vmax)
      span = ZoomSpan(axes, onselect, 'horizontal')

      onmove_callback is an optional callback that will be called on mouse move
      with the span range

    """

    def __init__(self, parent, axes, 
                        button=1,
                        minspan=None, 
                        useblit=False, 
                        rectprops=None, 
                        do_zoom_select_event=False, 
                        do_zoom_motion_event=False):
        """
        Create a span selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(vmin, vmax)

        If minspan is not None, ignore events smaller than minspan

        The span rect is drawn with rectprops; default
          rectprops = dict(facecolor='red', alpha=0.5)

        set the visible attribute to False if you want to turn off
        the functionality of the span selector


        """
        if rectprops is None:
            rectprops = dict(facecolor='yellow', alpha=0.2)

        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []

        self.rect = []
        self.background = None
        self.pressv = None
        self.axes_index = None

        self.rectprops = rectprops
        self.do_zoom_select_event = do_zoom_select_event
        self.do_zoom_motion_event = do_zoom_motion_event
        self.useblit = useblit
        self.minspan = minspan
        self.button = button

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prev = (0, 0)

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
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
        
        for axes in self.axes:
            trans = blended_transform_factory(axes.transData, axes.transAxes)
            self.rect.append(Rectangle( (0,0), 0, 1,
                                   transform=trans,
                                   visible=False,
                                   **self.rectprops ))

        if not self.useblit: 
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)


    def update_background(self, event):
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
        correct_button = event.button==self.button and not self.parent.shift_is_held
        return  event.inaxes not in self.axes or not self.visible or not correct_button

    def press(self, event):
        '''on button press event'''
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = False

        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i

        # remove the dynamic artist(s) from background bbox(s)
        for axes, rect in zip(self.axes, self.rect):
            if rect in axes.patches:
                axes.patches.remove(rect)
                self.canvas.draw()

        for rect in self.rect:
            rect.set_visible(self.visible)
            
        self.pressv = event.xdata
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = True

        for rect in self.rect:
            rect.set_visible(False)

        # left-click in place resets the x-axis 
        if event.xdata == self.pressv:
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds  
            xdel = self.parent.xscale_bump*(x1-x0)
            ydel = self.parent.yscale_bump*(y1-y0) / 1.1
            for axes in self.axes:
                if self.parent.reversex:
                    axes.set_xlim(x0+x1+xdel,x0-xdel)
                else:
                    axes.set_xlim(x0-xdel,x0+x1+xdel)
                # ylim is 0.1 at bottom to keep a 10:1 ratio between plot above
                # the zeroline and below, keeps zeroline from jittering. This
                # 10:1 ratio is set by set_ylim() method when the zeroline is 
                # at the top or bottom.
                axes.set_ylim(y0-ydel*0.1,y0+y1+ydel)
            self.canvas.draw()
            
            if self.do_zoom_select_event:
                self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True, iplot=self.axes_index)
            
            return

        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin>vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin
        if self.minspan is not None and span<self.minspan: return

        for axes in self.axes:
            if self.parent.reversex:
                axes.set_xlim((vmax, vmin))
            else:
                axes.set_xlim((vmin, vmax))
                #axes.set_xlim((self.pressv,event.xdata))
        self.canvas.draw()

        data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event, reversex=self.parent.reversex)
            self.parent.on_zoom_select(vmin, vmax, value, None, None, iplot=self.axes_index) 
            
        self.axes_index = None
        self.pressv = None
        
        return False

    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, rect in zip(self.axes, self.rect):
                axes.draw_artist(rect)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

    def onmove(self, event):
        '''on motion notify event'''
        if self.pressv is None or self.ignore(event): return
        x, y = event.xdata, event.ydata
        self.prev = x, y

        minv, maxv = x, self.pressv
        if minv>maxv: minv, maxv = maxv, minv
        for rect in self.rect:
            rect.set_x(minv)
            rect.set_width(maxv-minv)

        data_test = event.inaxes.lines!=[]

        if self.do_zoom_motion_event and data_test: 
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin
            value = self.parent.get_values(event, reversex=self.parent.reversex)
            self.parent.on_zoom_motion(vmin, vmax, value, None, None, iplot=self.axes_index) 

        self.update()
        return False




class CursorSpan:
    """
    Indicate two vertical reference lines along a matplotlib Axes

    Example usage:

      axes = subplot(111)
      axes.plot(x,y)

      def onselect(vmin, vmax):
          print(vmin, vmax)
      span = CursorSpan(axes, onselect)

      onmove_callback is an optional callback that will be called on mouse move
      with the span range

    """

    def __init__(self, parent, axes, 
                       button=3,
                       minspan=None, 
                       useblit=False, 
                       rectprops=None, 
                       do_refs_select_event=False, 
                       do_refs_motion_event=False):
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
        if rectprops is None:
            rectprops = dict(facecolor='none')

        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []

        self.rect = []
        self.background = None
        self.pressv = None
        self.axes_index = None

        self.rectprops = rectprops
        self.do_refs_select_event = do_refs_select_event
        self.do_refs_motion_event = do_refs_motion_event
        self.useblit = useblit
        self.minspan = minspan
        self.button = button

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prev  = (0,0)

        self.new_axes(axes)


    def set_span(self, xmin, xmax):
        
        x0, y0, x1, y1 = self.axes[0].dataLim.bounds

        self.visible = True
        if xmin < x0: xmin = x0
        if xmax > (x0+x1): xmax = x0+x1

        for rect in self.rect:
            rect.set_x(xmin)
            rect.set_width(xmax-xmin)
            
        self.canvas.draw()
        

    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
        
        for axes in self.axes:
            trans = blended_transform_factory(axes.transData, axes.transAxes)
            self.rect.append(Rectangle( (0,0), 0, 1,
                                   transform=trans,
                                   visible=False,
                                   **self.rectprops ))

        if not self.useblit: 
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)


    def update_background(self, event):
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
        correct_button = event.button==self.button or \
                         (event.button==self.parent.zoom_button and self.parent.shift_is_held)
        return  event.inaxes not in self.axes or not self.visible or not correct_button

    def press(self, event):
        '''on button press event'''
        self.visible = True
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_refs_motion_event:
            self.parent.do_motion_event = False

        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i

        # remove the dynamic artist(s) from background bbox(s)
        for axes, rect in zip(self.axes, self.rect):
            if rect in axes.patches:
                axes.patches.remove(rect)
                self.canvas.draw()
        
        for rect in self.rect:
            rect.set_visible(self.visible)
        self.pressv = event.xdata
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_refs_motion_event:
            self.parent.do_motion_event = True

        # left-click in place turns off display of cursors
        if event.xdata == self.pressv:
            self.visible = not self.visible
            for axes, rect in zip(self.axes, self.rect):
                rect.set_visible(self.visible)
                axes.add_patch(rect)
            self.canvas.draw()  
            
            if self.do_refs_select_event:
                data_test = event.inaxes.lines!=[]
                if self.do_refs_select_event and data_test:
                    vmin = self.pressv
                    rindex = get_data_index(event, xvalue=vmin, reversex=self.parent.reversex)
                else:
                    rindex = None
                    
                value = self.parent.get_values(event, reversex=self.parent.reversex)
                self.parent.on_refs_select(self.pressv, event.xdata, value, reset=True, iplot=self.axes_index, rindex=rindex)
            
            self.pressv = None  
            return
        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin>vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin
        # don't add reference span, if min span not achieved
        if self.minspan is not None and span<self.minspan: return

        for axes, rect in zip(self.axes, self.rect):
            rect.set_visible(True)
            axes.add_patch(rect)
            self.canvas.draw()

        data_test = event.inaxes.lines!=[]

        if self.do_refs_select_event and data_test:

            # update the reference lines data indices            
            imin = get_data_index(event, xvalue=vmin, reversex=self.parent.reversex)
            imax = get_data_index(event, xvalue=vmax, reversex=self.parent.reversex)
            if imin > imax: imin, imax = imax, imin
            self.parent.ref_locations = imin, imax
            
            # don't gather values if no onselect event
            value = self.parent.get_values(event, reversex=self.parent.reversex)
            self.parent.on_refs_select(vmin, vmax, value, reset=False, iplot=self.axes_index)
        
        self.axes_index = None
        self.pressv = None
        
        return False

    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, rect in zip(self.axes, self.rect):
                axes.draw_artist(rect)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

    def onmove(self, event):
        '''on motion notify event'''
        if self.pressv is None or self.ignore(event): return
        x, y = event.xdata, event.ydata
        self.prev = x, y

        minv, maxv = x, self.pressv
        if minv>maxv: minv, maxv = maxv, minv
        for rect in self.rect:
            rect.set_x(minv)
            rect.set_width(maxv-minv)

        data_test = event.inaxes.lines!=[]

        if self.do_refs_motion_event and data_test:
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin

            # update the reference lines data indices            
            imin = get_data_index(event, xvalue=vmin, reversex=self.parent.reversex)
            imax = get_data_index(event, xvalue=vmax, reversex=self.parent.reversex)
            if imin > imax: imin, imax = imax, imin
            self.parent.ref_locations = imin, imax

            # get data values at the current cursor location            
            value = self.parent.get_values(event, reversex=self.parent.reversex)            
            self.parent.on_refs_motion(vmin, vmax, value, iplot=self.axes_index) 

        self.update()
        return False




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
    def __init__(self, parent, axes, 
                             button=1,
                             drawtype='box',
                             minspanx=None, 
                             minspany=None, 
                             useblit=False,
                             lineprops=None, 
                             rectprops=None,
                             do_zoom_select_event=False, 
                             do_zoom_motion_event=False,
                             spancoords='data'):

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
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
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
        '''on button press event'''
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
        '''on button release event'''
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
                if self.parent.reversex:
                    axes.set_xlim(x0+x1+xdel,x0-xdel)
                else:
                    axes.set_xlim(x0-xdel,x0+x1+xdel)
                axes.set_ylim(y0-ydel,y0+y1+ydel)
            self.canvas.draw()
            
            if self.do_zoom_select_event:
                self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True, iplot=self.axes_index)
            
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
            if self.parent.reversex:
                axes.set_xlim((xmax,xmin))
            else:
                axes.set_xlim((xmin,xmax))
            axes.set_ylim((ymin,ymax))
        self.canvas.draw()


        data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event, reversex=self.parent.reversex)            
            self.parent.on_zoom_select(xmin, xmax, value, ymin, ymax, iplot=self.axes_index) # zeros are for consistency with box zoom

        self.axes_index = None
        self.eventpress   = None              # reset the variables to their
        self.eventrelease = None              #   inital values
        
        return False


    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
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
        '''on motion notify event if box/line is wanted'''
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

        data_test = event.inaxes.lines!=[]
        
        if self.do_zoom_motion_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event, reversex=self.parent.reversex)
            self.parent.on_zoom_motion(minx, maxx, value, miny, maxy, iplot=self.axes_index) # zeros are for consistency with box zoom
        
        self.update()
        return False

    def set_active(self, active):
        ''' Use to de/activate RectangleSelector with a boolean variable 'active' '''
        self.active = active

    def get_active(self):
        ''' to get status of active mode (boolean variable)'''
        return self.active



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
                       button=None,
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
        self.button = button

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
        '''return True if event should be ignored'''
        return  event.inaxes not in self.axes or event.button != self.button

    def press(self, event):
        '''on button press event'''
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
            self.parent.on_middle_press(event.x, event.y, self.axes_index, event.xdata, event.ydata)
        
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressxy is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = True

        xstr, ystr = self.pressxy
        xend = event.x #or self.prev[0]
        yend = event.y #or self.prev[1]

        if self.do_middle_select_event:
            self.parent.on_middle_select(xstr, ystr, xend, yend, self.axes_index)

        self.axes_index = None
        self.pressxy = None
        return False

    def onmove(self, event):
        '''on motion notify event'''
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
                menu.AppendSubMenu(subMenu, label)
            else:
                self.createMenuItem(self2, menu, *eachItem)
        return menu

    def createMenuItem(self, self2, menu, label, status, handler, kind=wx.ITEM_NORMAL):
        if not label:
            menu.AppendSeparator()
            return
        menuItem = menu.Append(-1, label, status, kind)
        self2.Bind(wx.EVT_MENU, handler, menuItem)




class DemoPlotPanel(PlotPanelScatter):
    """Plots several lines in distinct colors."""

    # Activate event messages
    _EVENT_DEBUG = True
    
    def __init__( self, parent, tab, **kwargs ):
        # initiate plotter
        PlotPanelScatter.__init__( self, parent, **kwargs )  
        self.tab    = tab
        self.top    = wx.GetApp().GetTopWindow()
        self.parent = parent
        self.count = 0

    def on_motion(self, xdata, ydata, value, bounds, iaxis):
        self.top.statusbar.SetStatusText( " Points = %i" % (xdata-0.5, ), 0)
        self.top.statusbar.SetStatusText( " Xval = %.2f" % (xdata, ), 1)

    def on_scroll(self, button, step, iaxis):
        self.set_vertical_scale(step)

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        delta  = xmax - xmin
        self.top.statusbar.SetStatusText(( " Point Range = %.2f to %.2f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText(( " dPoints = %i " % (delta, )), 2)

    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None, rindex=None):
        # Calculate area of span
        all_areas, all_rms = self.calculate_area()
        area = all_areas[0]
        rms  = all_rms[0]
        self.top.statusbar.SetStatusText(' Area = %1.5g  RMS = %1.5g' % (area,rms), 3)

        
    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        delta  = xmax - xmin
        self.top.statusbar.SetStatusText(( " Point Range = %.2f to %.2f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText(( " dPoints = %i " % (delta, )), 2)

        all_areas, all_rms = self.calculate_area()
        area = all_areas[0]
        rms  = all_rms[0]
        self.top.statusbar.SetStatusText(' Area = %1.5g  RMS = %1.5g' % (area,rms), 3)
        


class MyFrame(wx.Frame):
    def __init__(self, title="New Title Please", size=(350,200)):

        wx.Frame.__init__(self, None, title=title, pos=(150,150), size=size)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        util_CreateMenuBar(self)
        
        self.statusbar = self.CreateStatusBar(4, 0)

        self._create_fake_prefs()
        
        lines = self.make_data(points=90, lines=18)
        data = [[lines], [lines], [lines]]
        
        self.nb = wx.Notebook(self, -1, style=wx.BK_BOTTOM)
        
        panel1 = wx.Panel(self.nb, -1)
        self.view = DemoPlotPanel( panel1, self, 
                                      naxes=3,
                                      reversex=True,
                                      zoom='span', 
                                      reference=True,
                                      middle=True,
                                      unlink=False,
                                      do_zoom_select_event=True,
                                      do_zoom_motion_event=True,
                                      do_refs_select_event=True,
                                      do_refs_motion_event=True,
                                      do_middle_select_event=True,
                                      do_middle_motion_event=True,
                                      do_middle_press_event=False,
                                      do_scroll_event=True,
                                      xscale_bump=0.0,
                                      yscale_bump=0.05,
                                      data = None, 
                                      prefs=self.prefs, 
                                      )
    
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        panel1.SetSizer(sizer)
        self.view.Fit()    
    
        self.nb.AddPage(panel1, "One")
        
        self.view.set_data(data)
    
        self.view.set_color( (255,255,255) )


    def menuData(self):
        return [("&File", (
                    ("", "", ""),
                    ("&Quit",    "Quit the program",  self.on_close))),
                ("View", (
                    ("Zero Line", (
                        ("&Show",    "", self.on_zero_line_show,   wx.ITEM_CHECK),
                        ("", "", ""),
                        ("&Top",     "", self.on_zero_line_top,    wx.ITEM_RADIO),
                        ("&Middle",  "", self.on_zero_line_middle, wx.ITEM_RADIO),
                        ("&Bottom",  "", self.on_zero_line_bottom, wx.ITEM_RADIO))),
                    ("", "", ""),
                    ("X-Axis - Show", "", self.on_xaxis_show,       wx.ITEM_CHECK),
                    ("", "", ""),
                    ("Data Type - Summed", "", self.on_data_type_summed,    wx.ITEM_CHECK),
                    ("", "", ""),
                    ("&Placeholder",    "non-event",  self.on_placeholder))),
                ("Tests", (
                    ("Data 90 pts 18 lines",  "", self.on_test_90pts_18lines,  wx.ITEM_RADIO),
                    ("Data 50 pts 18 lines",  "", self.on_test_50pts_18lines,  wx.ITEM_RADIO),
                    ("Data 90 pts  8 lines",  "", self.on_test_90pts_8lines,   wx.ITEM_RADIO),
                    ("Data 50 pts  8 lines",  "", self.on_test_50pts_8lines,   wx.ITEM_RADIO),
                    ("", "", ""),
                    ("Show all Three", "", self.on_show_three, wx.ITEM_RADIO),
                    ("Show only One",  "", self.on_show_one,   wx.ITEM_RADIO),
                    ("", "", ""),
                    ("Placeholder",    "non-event",  self.on_placeholder)))]    


    def make_data(self, lines=8, points=50):

        if lines > 18: lines = 18
        if lines < 1:  lines = 1

        # create some data
        nlin = lines
    
        lines = np.mod(np.arange(nlin*points).reshape(nlin,points), points) - (float(points)/2.0)
        
        return lines


    def on_close(self, event):
        dlg = wx.MessageDialog(self, 
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()
            
    def on_zero_line_show(self, event):
        # won't need these with a real Prefs module
        self.prefs.zero_line_plot_show = not self.prefs.zero_line_plot_show
        # will need these
        self.view.update_axes()
        self.view.canvas.draw()

    def on_zero_line_top(self, event):
        # won't need these with a real Prefs module
        self.prefs.zero_line_plot_top = True
        self.prefs.zero_line_plot_middle = False
        self.prefs.zero_line_plot_bottom = False
        # will need these
        self.view.update_axes()
        self.view.canvas.draw()

    def on_zero_line_middle(self, event):
        # won't need these with a real Prefs module
        self.prefs.zero_line_plot_top = False
        self.prefs.zero_line_plot_middle = True
        self.prefs.zero_line_plot_bottom = False
        # will need these
        self.view.update_axes()
        self.view.canvas.draw()

    def on_zero_line_bottom(self, event):
        # won't need these with a real Prefs module
        self.prefs.zero_line_plot_top = False
        self.prefs.zero_line_plot_middle = False
        self.prefs.zero_line_plot_bottom = True
        # will need these
        self.view.update_axes()
        self.view.canvas.draw()

    def on_xaxis_show(self, event):
        # won't need these with a real Prefs module
        self.prefs.xaxis_show = not self.prefs.xaxis_show
        # will need these
        self.view.update_axes()
        self.view.canvas.draw()

    def on_data_type_summed(self, event):
        # will need these
        self.view.set_data_type_summed()
        self.view.update_plots()
        self.view.canvas.draw()

    def on_test_90pts_18lines(self, event):
        lines = self.make_data(lines=18, points=90)
        data = [[lines],[lines],[lines]]
        self.view.set_data(data)
        self.view.update(no_draw=False)

    def on_test_50pts_18lines(self, event):
        lines = self.make_data(lines=18, points=50)
        data = [[lines],[lines],[lines]]
        self.view.set_data(data)
        self.view.update(no_draw=False)

    def on_test_90pts_8lines(self, event):
        lines = self.make_data(lines=8, points=90)
        bob = lines.copy()
        dbob  = {'data':bob,'line_color_real':'green', 'markevery':[12, 17, 18, 19], 'markevery_color':'purple'}
        dbob2 = {'data':bob, 'markevery':[12, 17, 18, 19], 'markevery_color':'red'}
        data = [[lines],[dbob],[dbob2]]
        self.view.set_data(data)
        self.view.update(no_draw=False)

    def on_test_50pts_8lines(self, event):
        lines = self.make_data(lines=8, points=50)
        data = [[lines],[lines],[lines]]
        self.view.set_data(data)
        self.view.update(no_draw=False)

    def on_placeholder(self, event):
        print("Event handler for on_placeholder - not implemented")

    def on_show_one(self, event):
        self.view.change_naxes(1)

    def on_show_three(self, event):
        self.view.change_naxes(3)

    def _create_fake_prefs(self):
        self.prefs = fake_prefs()
    

#------------------------------------------------------------------------------

class fake_prefs(object):

    def __init__(self):

        self.foreground_color = "black"
        self.bgcolor = "#ffffff"
        self.zero_line_plot_show = False
        self.zero_line_plot_top = False
        self.zero_line_plot_middle = False
        self.zero_line_plot_bottom = True
        self.xaxis_show = False
        self.data_type_summed = False
        self.zero_line_plot_color = "goldenrod"
        self.zero_line_plot_style = "solid"
        self.line_color_real = "blue"
        self.line_color_imaginary = "red"
        self.line_color_magnitude = "purple"
        self.line_width      = 1.0     
        self.plot_view = "all"
#------------------------------------------------------------------------------

if __name__ == '__main__':

    app   = wx.App( 0 )
    frame = MyFrame( title='WxPython and Matplotlib', size=(600,600) )
    frame.Show()
    app.MainLoop()


"""
Expansion of matplotlib embed in wx example by John Bender and Edward 
Abraham, see http://www.scipy.org/Matplotlib_figure_in_a_wx_panel

This version allows the user to zoom in on the figure using either 
a span selector or a box selector. You can also set a persistent span
selector that acts as cursor references on top of whatever is plotted

ZoomSpan based on matplotlib.widgets.SpanSelector
CursorSpan based on matplotlib.widgets.SpanSelector
BoxZoom based on matplotlib.widgets.RectangleSelector

Brian J. Soher, Duke University, 20 October, 2010
"""

# Python modules

import math

# 3rd party modules
import matplotlib
import wx


# If we set the backend unconditionally, we sometimes get an undesirable
# message. 
if matplotlib.get_backend() != "WXAgg":
    matplotlib.use('WXAgg')

from matplotlib.transforms import blended_transform_factory
from matplotlib.patches    import Rectangle
from matplotlib.lines      import Line2D

# Our modules


    


class PlotPanel(wx.Panel):
    """
    The PlotPanel has a Figure and a Canvas and 'n' Axes. The user defines
    the number of axes on Init and this number cannot be changed thereafter.
    However, the user can select which of the n axes are actually displayed 
    in the Figure.
    
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
                               uses_collections=False,
                               xscale_bump=0.0,
                               yscale_bump=0.0,
                               props_zoom=None,
                               props_cursor=None,
                               update_plot=False,
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
        self.uses_collections = uses_collections        # has to be either Lines OR Collections, not both
        self.unlink = unlink
        self.xscale_bump = xscale_bump
        self.yscale_bump = yscale_bump

        self.update_plot = update_plot
        
        # Under GTK track self's size to avoid a continuous flow of size events.
        self._platform_is_gtk = ("__WXGTK__" in wx.PlatformInfo)
        self._current_size = (-1, -1)

        self.figure = Figure( None, dpi )
        #self.figure.subplots_adjust(top=0.98, right=0.98, hspace=0.08)

        # 'tight_layout' had a bug that shifted the right/top subplot locs really tight
        #   when I zoom and shifted back looser when I click in place to reset. Tried to
        #   bug hunt, but fell back on absolute subplots_adjust for now. bjs 8/2020
        #
        #self.figure.set_tight_layout({'pad': 0.1, 'h_pad': 0.1, 'w_pad': 0.1})
        #self.figure.set_tight_layout(True)

        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )

        # Create the N axes, add to figure, also keep a permanent reference to
        # each axes so they can be added or removed from the figure interactively.

        self.axes   = []                    # dynamic list for display
        for i in range(naxes):
            self.axes.append(self.figure.add_subplot(naxes,1,i+1))
        self.naxes = naxes   
        self.all_axes = list(self.axes)     # static list for long term interaction

        self.zoom = []
        self.refs = []
        self.middle = []
        
        self.set_color( color )
        self._set_size()
        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._on_idle)
        self.Bind(wx.EVT_SIZE, self._on_size)

        # ensure that default properties exist
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
                                           do_middle_select_event=do_middle_select_event,
                                           do_middle_motion_event=do_middle_motion_event,
                                           do_middle_press_event=do_middle_press_event)
            else:
                for axes in self.axes:
                    self.middle.append( MiddleEvents(self, [axes],
                                          button=middle_button,
                                          do_middle_select_event=do_middle_select_event,
                                          do_middle_motion_event=do_middle_motion_event))

        self.do_motion_event = True  
        self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)
        
        self.do_scroll_event = do_scroll_event
        if self.do_scroll_event:
            self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)




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
        if self.update_plot:
            self.on_update_plot()

    def _on_move(self, event):
        """
        Internal method organizes data sent to an external user defined event
        handler for motion events. We gather values from either collection or
        line plots, determine which axis we are in, then call on_motion()

        """
        if event.inaxes == None or not self.do_motion_event: return
        
        bounds = event.inaxes.dataLim.bounds
        value = self.get_values(event)
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i

        self.on_motion(event.xdata, event.ydata, value, bounds, iaxis)
        
    def _on_scroll(self, event):
        """
        Internal method organizes data sent to an external user defined event
        handler for scroll events. We determine which axis we are in, then
        call on_scroll()

        """
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i
                
        self.on_scroll(event.button, event.step, iaxis)        


    def update_plot(self):
        # placeholder for user derived method
        pass

    def get_values(self, event):
        """
        Included in plot_panel object so users can overwrite it if necessary.
        Default functionality: Determine which axes the mouse is in, return a
        list of data values at the x location of the cursor.
        
        This is complicated by the fact that some plot_panels put their plots
        into the lines attribute, while others use the collections attribute.
        This is solved by designating this distinction using a flag at the 
        main level. Caveat, plot_panels can have only one or the other types 
        of plots, line or collection based.
        
        """
        value = []
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds

        if not all([math.isfinite(x) for x in [x0, y0, x1, y1]]):    # typical when plot is empty
            return [0.0,]
    
        # sporadic event calls yield 'event.inaxes == None' so we test for that here
        if event.inaxes:
            if self.uses_collections:
                
                # Expects LineCollection - Figures out the offset of the cursor and
                # returns the path value for the x-position of that line plot.
                # Note. Not guaranteed to work for any other types of collections
                if event.inaxes.collections:
                    npts = len(event.inaxes.collections[0]._paths[0].vertices[:,1])
                    indx = int(round((npts-1) * (event.xdata-x0)/x1))
                    if self.reversex:   indx = npts - indx - 1
                    if indx > (npts-1): indx = npts - 1
                    if indx < 0:        indx = 0
                    for i,path in enumerate(event.inaxes.collections[0]._paths):
                        offset = event.inaxes.collections[0]._uniform_offsets[i,1]
                        dat = path.vertices[:,1]
                        if indx < len(dat):
                            value.append(dat[indx]-offset)
            else:
                if event.inaxes.lines:
                    data = event.inaxes.lines[0].get_ydata()
                    if len(data.shape)==0:
                        value = [data]
                    else:
                        npts = len(event.inaxes.lines[0].get_ydata())
                        indx = int(round((npts-1) * (event.xdata-x0)/x1))
                        if self.reversex:   indx = npts - indx - 1
                        if indx > (npts-1): indx = npts - 1
                        if indx < 0:        indx = 0
                        for line in event.inaxes.lines:
                            dat = line.get_ydata()
                            if indx < len(dat):
                                value.append(dat[indx])
        if value == []: 
            value = [0.0,]
            
        return value


    #=======================================================
    #
    #           User Accessible Functions  
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

    def on_refs_select(self, xmin, xmax, val, iplot=None):
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
        self._dprint('on_middle_move, X(cur,prev)='+str(xstr)+','+str(xprev)+'  Y(cur,prev)='+str(ystr)+','+str(yprev)+'  Index = '+str(indx))

    def on_middle_press(self, xloc, yloc, indx, bounds=None, xdata=None, ydata=None):
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
        return  event.inaxes not in self.axes or not self.visible or event.button != self.button

    def press(self, event):
        '''on button press event'''
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = False

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
            if all([math.isfinite(x) for x in [x0, y0, x1, y1]]):  # typical when plot is empty
                xdel = self.parent.xscale_bump * (x1 - x0)
                ydel = self.parent.yscale_bump * (y1 - y0)
                for axes in self.axes:
                    if self.parent.reversex:
                        axes.set_xlim(x0+x1+xdel,x0-xdel)
                    else:
                        axes.set_xlim(x0-xdel,x0+x1+xdel)
                    axes.set_ylim(y0-ydel,y0+y1+ydel)
                self.canvas.draw()

                if self.do_zoom_select_event:
                    self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True)
            
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

        self.canvas.draw()

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes is not None:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event)
            self.parent.on_zoom_select(vmin, vmax, value, None, None) 
            
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

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_motion_event and data_test: 
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin
            value = self.parent.get_values(event)
            self.parent.on_zoom_motion(vmin, vmax, value, None, None) 

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
        if all([math.isfinite(x) for x in [x0, y0, x1, y1]]):  # typical when plot is empty
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
        return  event.inaxes not in self.axes or not self.visible or event.button != self.button

    def press(self, event):
        '''on button press event'''
        self.visible = True
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_refs_motion_event:
            self.parent.do_motion_event = False

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

        # left-click in place resets the x-axis 
        if event.xdata == self.pressv:
            self.visible = not self.visible
            for axes, rect in zip(self.axes, self.rect):
                rect.set_visible(self.visible)
                axes.add_patch(rect)
            self.canvas.draw()  
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

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_refs_select_event and data_test:
            # don't gather values if no onselect event
            value = self.parent.get_values(event)
            self.parent.on_refs_select(vmin, vmax, value)
        
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

        # sporadic event calls yield event.inaxes == None so we test here
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_refs_motion_event and data_test:
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin
            value = self.parent.get_values(event)            
            self.parent.on_refs_motion(vmin, vmax, value) 

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


        # make the drawn box/line visible get the click-coordinates, button, ...
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
            if all([math.isfinite(x) for x in [x0,y0,x1,y1]]):  # typical when plot is empty
                xdel = self.parent.xscale_bump * (x1 - x0)
                ydel = self.parent.yscale_bump * (y1 - y0)
                for axes in self.axes:
                    if self.parent.reversex:
                        axes.set_xlim(x0+x1+xdel,x0-xdel)
                    else:
                        axes.set_xlim(x0-xdel,x0+x1+xdel)
                    axes.set_ylim(y0-ydel,y0+y1+ydel)
                self.canvas.draw()

                if self.do_zoom_select_event:
                    self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True)
            
            return

        self.eventrelease = event   # release coordinates, button, ...

        if self.spancoords=='data':
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
            # calculate dimensions of box or line get values in the right order
        elif self.spancoords=='pixels':
            xmin, ymin = self.eventpress.x, self.eventpress.y
            xmax, ymax = self.eventrelease.x, self.eventrelease.y
        else:
            raise ValueError('spancoords must be "data" or "pixels"')

        # assure that min<max values and x and y values are not equal
        if xmin>xmax: xmin, xmax = xmax, xmin
        if ymin>ymax: ymin, ymax = ymax, ymin
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

        # sporadic event calls yield event.inaxes == None so we test here
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            value = self.parent.get_values(event)
            self.parent.on_zoom_select(xmin, xmax, value, ymin, ymax) # zeros are for consistency with box zoom

        self.eventpress   = None              # reset variables to inital values
        self.eventrelease = None
        
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

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]
        
        if self.do_zoom_motion_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event)
            self.parent.on_zoom_motion(minx, maxx, value, miny, maxy) # zeros are for consistency with box zoom
        
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
          print(vmin, vmax)
      middle = MiddleEvents(axes, onselect)

      onmove_callback is an optional callback that will be called on mouse move
      with the span range

    """

    def __init__(self, parent, axes,
                       button=2,
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
        self.button = button

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
        '''return True if event should be ignored'''
        return  event.inaxes not in self.axes or event.button !=self.button

    def press(self, event):
        '''on button press event'''
        if self.ignore(event): return
        self.buttonDown = True
        
        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = False

        bounds = event.inaxes.dataLim.bounds

        self.pressxy = event.x, event.y
        self.prevxy  = event.x, event.y

        if self.do_middle_press_event:
            self.parent.on_middle_press(event.x, event.y, self.axes_index, bounds=bounds, xdata=event.xdata, ydata=event.ydata)
        
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
        xend = event.x
        yend = event.y

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

if __name__ == '__main__':

    import numpy as np

    class DemoPlotPanel(PlotPanel):
        """Plots several lines in distinct colors."""

        # Activate event messages
        _EVENT_DEBUG = True

        def __init__( self, parent, **kwargs ):
            # initiate plotter
            PlotPanel.__init__( self, parent, **kwargs )  
            self.parent = parent


    app   = wx.App( False )
    frame = wx.Frame( None, wx.ID_ANY, 'WxPython and Matplotlib - PlotPanel', size=(800,800) )
    
    nb = wx.Notebook(frame, -1, style=wx.BK_BOTTOM)
    
    panel1 = wx.Panel(nb, -1)
    view = DemoPlotPanel( panel1, naxes=2,
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
                                  do_scroll_event=True )

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
    panel1.SetSizer(sizer)
    view.Fit()    

    nb.AddPage(panel1, "One")
    
    t = np.arange(0.0, 2.0, 0.01)
    s = 1.0*np.sin(2*np.pi*t)
    c = 1.0*np.cos(2*np.pi*t)

    view.set_color( (255,255,255) )

    line1 = view.axes[0].plot(t, s, 'b', linewidth=1.0)
    line2 = view.axes[1].plot(t, c, 'g', linewidth=3.0)
    view.axes[0].set_xlabel('time (s)')
    for axes in view.axes:
        axes.set_ylabel('voltage (mV)')
        axes.grid(True)
    view.axes[0].set_title('About as simple as it gets, folks')
    
    view.canvas.figure.subplots_adjust(left=0.20,
                                         right=0.95,
                                         bottom=0.15,
                                         top=0.95,
                                         wspace=0.0,
                                         hspace=0.01)
    view.canvas.draw()
    frame.Show()
    app.MainLoop()




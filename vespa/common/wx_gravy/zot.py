
import matplotlib
from matplotlib.patches import Circle
from matplotlib.widgets import AxesWidget


class Zot(AxesWidget):
    """
    A red star that can be plotted anywhere on top of an image.

      *horizOn*
        Controls the visibility of the horizontal line

      *vertOn*
        Controls the visibility of the horizontal line

    and the visibility of the cursor itself with the *visible* attribute
    """
    def __init__(self, ax, radius=10, useblit=False, **lineprops):
        """
        Add a cursor to *ax*.  If ``useblit=True``, use the backend-
        dependent blitting features for faster updates (GTKAgg
        only for now).  *lineprops* is a dictionary of line properties.

        .. plot :: mpl_examples/widgets/cursor.py
        """
        # TODO: Is the GTKAgg limitation still true?
        AxesWidget.__init__(self, ax)

        self.visible = True
        self.radius  = radius
        self.useblit = useblit and self.canvas.supports_blit

        if self.useblit:
            lineprops['animated'] = True
        self.zot = Circle((0,0), radius)

        self.background = None
        self.needclear = False


    def clearzot(self, event):
        """clear the cursor"""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.zot.set_visible(False)


    def setzot(self, xloc, yloc):
        """draw the zot somewhere if visible"""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.zot.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        self.needclear = True
        if not self.visible:
            return
        self.zot.set_xydata((xloc, yloc))
        self.zot.set_visible(self.visible)

        self._update()


    def _update(self):

        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.zot)
            self.canvas.blit(self.ax.bbox)
        else:

            self.canvas.draw_idle()

        return False
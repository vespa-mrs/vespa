# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.constants as constants

def is_oc(tab_or_transform):
    """
    Given a tab or a transform, returns True if the tab's
    transform (or the transform itself) is an optimal control
    transform, False otherwise.

    This isn't particularly tricky, it just encapsulates a question we need
    to ask in a few different places.
    """
    type_ = None

    # All of our tabs inherit from wx.Panel
    if isinstance(tab_or_transform, wx.Panel):
        # It's a tab
        tab = tab_or_transform
        if hasattr(tab, "transform") and tab.transform:
            type_ = tab.transform.type
    else:
        # It's a transform
        type_ = tab_or_transform.type

    return (type_ == constants.TransformationType.OCN)

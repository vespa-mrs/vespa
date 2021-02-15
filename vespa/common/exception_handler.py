# Python modules

import platform
import sys
import os
import urllib.request, urllib.parse, urllib.error
import webbrowser
import traceback as traceback_module

# 3rd party modules
import wx

# Our modules
import vespa.common.auto_gui.exception_report as gui_exception_report
import vespa.common.wx_gravy.util as wx_util
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time

"""A specialized module that contains our exception hook and the functions
that support it. An exception hook allows us to execute arbitrary code when
an unhandled exception occurs. Our hook is a superset of Python's standard
exception handling.

For more info on exception hooks, see the Python standard lib doc for 
sys.excepthook.
"""

# _EMAIL_ADDRESS is the address to which bug reports should be sent.
_EMAIL_ADDRESS = "vespa.bugs@gmail.com"

# _DEFAULT_VERSION is the string reported when we can't get version info for
# a module in the debug info.
_DEFAULT_VERSION = "???"

_DIALOG_TEXT = "Sorry, an unexpected problem occurred.\n\n"                  \
"To help us fix the problem, please email a report to us.\n\n"              \
"If you want to send the report via Web mail (like GMail or Yahoo! Mail), " \
"copy and paste the details into an email and send them to %s.\n\n"         \
"If you prefer to use a local email program, you can check the \"open a "   \
"new email\" checkbox to have %s create the email for you.\n\n"             \
"Once you've sent the message, we recommend that you restart %s."

# The dialog width is pretty arbitrary. 450ish looks nice.
_APPROX_DIALOG_WIDTH = 450

_HEADER = ("-" * 25) + " %s " + ("-" * 25)


class DialogExceptionReport(gui_exception_report.MyDialog):
    """Displays a dialog that informs them that something went awry, displays
    the exception details and offers an opportunity to email those details to
    us. Callers should pass the app name (e.g. "Simulation") and a string
    containing everything that the user should send in the problem report
    email (traceback + debug info).
    """

    def __init__(self, app_name, exception_details):
        gui_exception_report.MyDialog.__init__(self, wx.GetApp().GetTopWindow())

        self.app_name = app_name

        # Populate the message and explictly wrap the text to the width of
        # the dialog (leaving room for the borders). OS X wraps automatically,
        # Win & Gnome don't.
        label_text = _DIALOG_TEXT % (_EMAIL_ADDRESS, app_name, app_name)
        self.LabelMessage.SetLabel(label_text)
        self.LabelMessage.Wrap(_APPROX_DIALOG_WIDTH - 20)

        self.SetSize( (_APPROX_DIALOG_WIDTH, -1) )

        # Populate the textbox with the exception report. We use a textbox
        # so people can scroll and also select all/copy manually if they
        # like. On my version of OSX/wx, horizontal scrolling is disabled
        # on this readonly textbox which is unfortunate.
        textbox = self.TextDetails
        textbox.SetValue(exception_details)
        # Moving the insertion point to 0 ensures that the top of the text
        # is visible. (Works under OSX & Windows, not Ubuntu/Gnome.)
        textbox.SetInsertionPoint(0)

        # Set textbox font to monospace
        wx_util.set_font_to_monospace(textbox)

        # The dialog requires an explicit poke to account for the fact that
        # the message label has changed size. Under Windows, this call has
        # to come after all of the above for it to work properly.
        self.Fit()

        self.SetTitle("%s Exception Report" % app_name)

        self.Layout()

        self.Center()


    def on_copy(self, event):
        exception_report = self.TextDetails.GetValue()
        
        wx_util.copy_to_clipboard(exception_report)

        if self.CheckboxOpenEmail.IsChecked():
            _create_exception_email(self.app_name, exception_report)



def exception_hook(exception_type, value, traceback):
    """Our custom exception hook. See the Python standard library doc for
    sys.excepthook.
    """
    # We always invoke Python's standard exception handling first.
    sys.__excepthook__(exception_type, value, traceback)
    
    app_name = wx.GetApp().GetAppName()

    # Turn the traceback into a string.
    traceback = "\n".join(traceback_module.format_tb(traceback))
    traceback += "\n%s: %s" % (exception_type.__name__, value)
    
    # Next step is to log it to disk. I don't have any particular plans for
    # this file. I just write it to be thorough.
    filename = "vespa_exception.%s.txt" % util_time.filename_timestamp()
    filename = os.path.join(util_misc.get_data_dir(), "logs", filename)

    f = open(filename, "w")
    now = util_time.now(util_time.ISO_TIMESTAMP_FORMAT)
    s = "%s exception logged at %s (local time)\n" % (app_name, now)
                                
    s += ("-" * 60) + "\n"
    f.write(s)
    f.write(traceback)

    f.close()

    # Last step is to display the exception dialog to the user.
    lines = [_HEADER % (" " + app_name + " Exception ")]
    lines.append(_HEADER % " TRACEBACK ")
    lines.append(traceback)
    lines += ["*%s\t%s" % line for line in _debug_report()]
    s = "\n".join(lines) + "\n"

    dialog = DialogExceptionReport(app_name, s)
    dialog.ShowModal()
    dialog.Destroy()


def _create_exception_email(app_name, exception_report):
    """Constructs a mailto: URL that will (hopefully) open the user's default
    mail client.
    """
    # Adding a timestamp to the subject will help us keep error reports
    # separate in our inboxes.
    timestamp = util_time.now("%d %B %Y, %H:%M:%S")
    subject = "%s Exception Report at %s" % (app_name, timestamp)

    # Don't make the body text too long. Remember that the addressee, subject
    # and body all have to fit into a URL. The HTTP spec defines no limit on
    # the length of a URL, but the application handling the URL certainly has
    # a practical limit. We can't predict what that limit will be.
    # FWIW, IE8 refuses URLs much longer than 2000 characters. I suspect the
    # lowest common denominator mail app is a lot shorter. The body text
    # below produces a URL that approaches 256 characters.
    # ref: http://support.microsoft.com/kb/q208427/
    # ref: http://stackoverflow.com/questions/417142/what-is-the-maximum-length-of-an-url
    body = "Just hit paste and send. Optionally, you can also tell us what " \
           "you were doing when the error occurred."
    url = "mailto:%s?&subject=%s&body=%s" % (urllib.parse.quote(_EMAIL_ADDRESS),
                                             urllib.parse.quote(subject),
                                             urllib.parse.quote(body))

    webbrowser.open(url)


def _debug_report(reports=["platform", "python", "version"]):
    """Builds a report useful for troubleshooting. Contains info about
    the local machine.
    """
    lines = [ ]
    
    for report in reports:
        lines.append( (_HEADER % report.upper(), "") )

        if report == "platform":
            lines += _get_platform_info()
        elif report == "python":
            lines += _get_python_info()
        elif report == "version":
            lines += _get_version_info()

    return lines



def _get_module_version(module):
    """Given a module object (e.g. sys.modules[0]), returns the module's
    version string if it's available under one of several common attributes.
    """
    version = _DEFAULT_VERSION
    for attribute_name in ("__version__", "__VERSION__", "VERSION",
                           "version"):
        if hasattr(module, attribute_name):
            version = getattr(module, attribute_name)
            if callable(version):
                try:
                    version = version()
                except TypeError:
                    # module.version() requires params; there's not much we
                    # can do about that.
                    version = _DEFAULT_VERSION
            break

    return str(version)


def _get_platform_info():
    """Returns details about the platform (OS, patch level, etc.) in a
    list of 2-tuples of (name, value).
    """
    lines = [ ]
    lines.append( ("platform.platform", platform.platform()) )
    lines.append( ("platform.release", platform.release()) )
    lines.append( ("platform.uname", platform.uname()) )

    platform_name = sys.platform.lower()

    if platform_name.startswith("win"):
        # Under Windows, the platform name will be win32 or win64 I think.
        # AFAIK the call sys.getwindowsversion() should always be available
        # if the platform is Windows but I'd rather be on the safe side.
        if hasattr(sys, "getwindowsversion"):
            lines.append( ("sys.getwindowsversion()", str(sys.getwindowsversion())))

    if "darwin" in platform_name:
        lines.append( ("platform.mac_ver()", str(platform.mac_ver())) )

    if "linux" in platform_name:
        # There are two calls to get Linux distro information. The first
        # is preferred but not present in Python 2.5. The latter is
        # deprecated.
        if hasattr(platform, "linux_distribution"):
            lines.append( ("platform.linux_distribution()",
                            str(platform.linux_distribution())) )
        elif hasattr(platform, "dist"):
            lines.append( ("platform.dist()", str(platform.dist())) )

    return lines


def _get_version_info():
    """Returns version info for Vespa and SQLite. Also attempts to return
    version info for all loaded Python modules.

    Info is returned in a list of 2-tuples of (name, value).
    """
    import sqlite3
    lines = [ ]

    lines.append( ("vespa version", util_misc.get_vespa_version()) )
    lines.append( ("sqlite3.sqlite_version", sqlite3.sqlite_version) )

    # Here I build a list of the currently loaded modules. The bulk of this
    # code is here to filter out noise. At any given moment there's lots of
    # modules loaded that we (mostly) don't care about. For instance, try
    # this to see how many modules are loaded just as a result of invoking
    # Python:
    #    python -c "import sys;  print sys.modules"
    #

    # Some don't have a __file__ attr
    modules = [module for module in list(sys.modules.values())
                                 if hasattr(module, "__file__")]

    # I figure out where most of the standard library modules are by selecting
    # a module (I picked the re module) and getting that module's path. I
    # then assume that anything with the same path is part of the Python
    # standard lib. This doesn't exclude all of the modules I'd like it to,
    # but it successfully removes many while erring on the side of caution.
    standard_library_path = ""
    for module in modules:
        head, tail = os.path.split(module.__file__)
        if tail == "re.pyc":
            standard_library_path = head

    # Remove Python modules
    f = lambda fq_filename: os.path.split(fq_filename)[0] != standard_library_path
    modules = [module for module in modules if f(module.__file__)]

    # What's left should be mostly modules that are not part of the standard
    # library. I build a list of 2-tuples consisting of (module filename,
    # module version).
    modules = [ (module.__file__, _get_module_version(module)) for module
                                                               in  modules]

    modules = sorted(modules)

    return lines + modules


def _get_python_info():
    """Returns info about Python (version, path to executable, etc.) in a
    list of 2-tuples of (name, value)."""
    lines = [ ]

    lines.append( ("bitiness", "%d-bit" % util_misc.get_bit_mode()) )
    lines.append( ("sys.version", sys.version) )
    lines.append( ("sys.executable", sys.executable) )
    lines.append( ("sys.path", sys.path) )

    return lines



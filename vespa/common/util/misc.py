"""Miscellaneous utilities"""

# Python modules

import os
import uuid as uuid_module
import sys
import re

try:
    import importlib.metadata 
    version_method = 'importlib'
except:
    # 3rd party imports
    import pkg_resources
    version_method = 'pkg_resources'

# Beware -- lots of other Vespa modules import this one. If you import
# other Vespa modules here you will almost certainly create circular imports!

# The _EOL_REGEX matches all possible line endings: \r\n (Windows),
# \r (old Macs) and \n (Unix).
# This is used by normalize_newlines().
_EOL_REGEX = re.compile(r"(?:\r\n)|\r|\n")

# _ISOTOPE_REGEX matches isotope strings, e.g. 1H, 31P, P31, C13, 7Li, or 29si.
# Not case-sensitive. Assumes all whitespace has been removed.
# This is used by normalize_isotope_name().
_ISOTOPE_REGEX = re.compile("""
    (\d{1,2})      # Match & capture 1 or 2 digits
    ([A-Z]{1,2})   # Match & capture 1 or 2 letters
    $              # String must end here
    |              # ...or...
    ([A-Z]{1,2})   # Match & capture 1 or 2 letters
    (\d{1,2})      # Match & capture 1 or 2 digits
    $              # String must end here
                           """, re.VERBOSE|re.IGNORECASE)

# _ALL_DIGITS is used by normalize_isotope_name().
_ALL_DIGITS = [str(i) for i in range(10)]

#############    Classes start here

class ClassProperty(object):
    """Creates a read-only property associated with a class (rather than
    a class instance). If Python had a @classproperty decorator, this is
    how it would work except that it's not possible to create a set-able
    property with this trick.

    I stole this code directly from Michael Foord:
    http://mail.python.org/pipermail/python-ideas/2011-January/008958.html

    To instantiate, pass a function that accepts one param (which will be
    the associated class). For example:

        def get_the_answer(klass):
            # The klass param is unused here
            return 42

        class DeepThought(object):
            ANSWER = ClassProperty(get_the_answer)

        print DeepThought.ANSWER
        42

    A shorter example, using lambda:
        class DeepThought(object):
            ANSWER = ClassProperty(lambda klass: 42)

        print DeepThought.ANSWER
        42
    """
    def __init__(self, function):
        self._function = function

    def __get__(self, instance, klass):
        return self._function(klass)


class StaticProperty(object):
    """Identical to ClassProperty (q.v.), except that the function passed
    to the constructor must take no arguments. For example:

        def get_the_answer():
            return 42

        class DeepThought(object):
            ANSWER = StaticProperty(get_the_answer)

        print DeepThought.ANSWER
        42


    A shorter example, using lambda:
        class DeepThought(object):
            ANSWER = StaticProperty(lambda: 42)

        print DeepThought.ANSWER
        42
    """
    def __init__(self, function):
        self._function = function

    def __get__(self, instance, klass):
        return self._function()


class WindowsSpecialFolderIds(object):
    """Constants for some Windows special folders. These are useful for
    passing to our function get_windows_special_folder_path().

    For details on special folders, see here:
    http://msdn.microsoft.com/en-us/library/bb762494%28v=vs.85%29.aspx
    """
    # If you want to use a value that's not defined here, just add it. CSIDL
    # values are defined in shlobj.h. Here's one source for that file:
    # http://source.winehq.org/source/include/shlobj.h
    CSIDL_DESKTOP       = 0x0000
    CSIDL_PERSONAL      = 0x0005
    CSIDL_MYDOCUMENTS   = CSIDL_PERSONAL
    CSIDL_APPDATA       = 0x001a
    CSIDL_LOCAL_APPDATA = 0x001c


#############    Functions start here

def get_bit_mode():
    """Returns 32 or 64 (as an int) to indicate whether Python is compiled
    as a 32- or 64-bit app. Note that this doesn't tell you whether or not
    the operating system is a 64-bit OS, it only informs about Python.
    """
    # It's trickier than you might think to get this information. Python's
    # platform.architecture() can get confused under OS X, and the
    # once-preferred alternative of testing sys.maxint doesn't work under
    # Win64.
    # The solution below is blessed by the wisdom of stackoverflow.com.
    # References --
    #     http://stackoverflow.com/questions/1405913/how-do-i-determine-if-my-python-shell-is-executing-in-32bit-or-64bit-mode
    #     http://stackoverflow.com/questions/3411079/why-does-the-python-2-7-amd-64-installer-seem-to-run-python-in-32-bit-mode
    #     http://mail.python.org/pipermail/python-list/2010-October/1258275.html
    #     http://groups.google.com/group/comp.lang.python/msg/5e363fcd9131dec4
    if hasattr(sys, "maxsize"):
        # This works under Python >= 2.6
        return 64 if (sys.maxsize > 2**32) else 32
    else:
        # This works under Python 2.5.
        import struct
        return 8 * struct.calcsize("P")


def get_data_dir():
    r"""Returns the path to the user's data folder which is where we store our
    config files and database, among other things. The path depends on the user
    name and operating system.

    This is the same as the path returned by
    wx.StandardPaths.Get().GetUserDataDir(), but this function isn't dependent
    on wx. The wx GetUserDataDir() function returns a path that contains the
    value returned by wx.GetAppName(). This call always uses "Vespa".

    You can expect data dir names something like this:
      * Windows: C:\Documents and Settings\%USERNAME%\Application Data\Vespa
      * Linux:   ~/.Vespa
      * OS X:    ~/Library/Application Support/Vespa
    """
    APP_NAME = "Vespa"

    data_dir = None

    home = os.path.expanduser("~")

    platform = get_platform()

    if platform == "linux":
        data_dir = os.path.join(home, "." + APP_NAME)

        if not os.path.exists(data_dir):
            # Up until now (March 2011), we've been putting our data files in
            # ~/.Vespa. To be better citizens, we'll now begin to respect the
            # new-ish XDG standard described on freedesktop.org. So on machine
            # where ~/.Vespa doesn't yet exist (i.e. new installs), we'll try
            # the location that the XDG standard recommends.
            data_dir = os.getenv("XDG_DATA_HOME", "")
            if os.path.exists(data_dir):
                data_dir  = os.path.join(data_dir, "." + APP_NAME)
            else:
                # Oh well, at least we tried.
                data_dir = os.path.join(home, "." + APP_NAME)
    elif platform == "osx":
        # Apple's doc says, "The preferred location for nearly all support
        # files is...the current user's ~/Library/Application Support
        # directory. Within the Application Support directory, you should
        # always place support files in a custom subdirectory named for your
        # application or company."
        # ref: http://developer.apple.com/library/mac/#documentation/MacOSX/Conceptual/BPFileSystem/Articles/WhereToPutFiles.html
        data_dir = os.path.join(home, "Library/Application Support", APP_NAME)
    elif platform == "windows":
        # Note that we use the local (non-roaming) application data directory.
        data_dir = get_windows_special_folder_path(WindowsSpecialFolderIds.CSIDL_LOCAL_APPDATA)
        data_dir = os.path.join(data_dir, APP_NAME)

    return data_dir


def get_documents_dir():
    """Returns the path to the user's "My Documents" or similarly-named
    folder. If it can't find one, it returns the user's home folder.

    When the Documents directory exists, this function should return the
    same value as the wx function GetDocumentsDir(). This function, however,
    isn't dependent on wx.
    """
    documents_dir = None

    home = os.path.expanduser("~")

    platform = get_platform()

    if platform == "linux":
        documents_dir = get_linux_special_folder_path("XDG_DOCUMENTS_DIR")
    elif platform == "osx":
        documents_dir = os.path.join(home, "Documents")
    elif platform == "windows":
        documents_dir = get_windows_special_folder_path(WindowsSpecialFolderIds.CSIDL_MYDOCUMENTS)

    if not documents_dir:
        documents_dir = home

    return documents_dir


def get_linux_special_folder_path(folder_name):
    """Given an XDG_XXX string (e.g. "XDG_DOCUMENTS_DIR"), attempts to read
    the path from user-dirs.dirs. If it can't find the file or the path, or
    the path is specified but doesn't exist, the function returns "".
    """
    # This code is based on the descriptions here:
    # http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
    # http://www.freedesktop.org/wiki/Software/xdg-user-dirs
    path = ""

    home = os.path.expanduser("~")

    xdg_config_home = os.getenv("XDG_CONFIG_HOME", "")

    if not xdg_config_home:
        xdg_config_home = os.path.join(home, ".config")

    dirs_filename = os.path.join(xdg_config_home, "user-dirs.dirs")

    if os.path.exists(dirs_filename):
        f = open(dirs_filename)

        for line in f:
            line = line.strip()

            if line.startswith(folder_name):
                # Line looks like this:
                # XDG_DOCUMENTS_DIR="$HOME/Documents"
                i = line.find("=")
                i = len(line) - i
                line = line[i:]
                line = line.strip('"')

                path = line.replace("$HOME", home)

                break
        f.close()

    if not os.path.exists(path):
        path = ""

    return path


def get_platform():
    """Returns the current platform/operating system as one of "osx", "linux"
    or "windows". The function is intentionally coarse-grained and doesn't
    differentiate between different versions of these operating systems.

    If the platform can't be determined, something is probably seriously
    wrong and this function returns None.
    """
    platform = sys.platform.lower()

    if "linux" in platform:
        simple_platform = "linux"
    # Beware! If you simply test for ("win" in platform) you will get a
    # surprise when "darwin" is identified as Windows.
    elif platform.startswith("win"):
        simple_platform = "windows"
    elif "darwin" in platform:
        simple_platform = "osx"
    else:
        simple_platform = None

    return simple_platform


def get_vespa_install_directory():
    """Returns the fully-qualified name of the directory in which Vespa has
    been installed.

    Under OS X this will return something like --
    /Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/vespa

    Under Linux this will return something like --
    /usr/local/lib/python2.6/dist-packages/vespa

    Under Windows this will return something like --
    c:\\Python25\\lib\\site-packages\\vespa
    """
    # We figure out the path by picking a file at a known location in the
    # vespa hierarchy, getting that file's path and inferring the rest.
    path = ""

    if hasattr(sys, "frozen") and sys.frozen:
        # Looks like this is being used in a PyInststaller bundle. This never
        # happens in Vespa, only when this module is used in the Siemens
        # realtime project.
#        path = os.environ['_MEIPASS2']
        path = sys._MEIPASS
    else:
        # This is the normal case in Vespa.

        # Strange but true -- we can 'import vespa', although it's meaningless
        # to do so in a practical sense. However, we can then use vespa.__file__
        # to see where that module resides.

        import vespa
        path = os.path.abspath(vespa.__file__)
        # path is now a fully-qualified pathname with "__init__.pyc"
        # dangling from the end.
        # Trim the filename
        path = os.path.dirname(path)

    return path


def get_vespa_version():
    """Returns Vespa's version as a string, e.g. "1.0.1"."""
    path = os.path.join(get_vespa_install_directory(), '..', 'VERSION')

    if os.path.exists(path):
        # The VERSION file isn't distributed as part of Vespa, so if I can see it, this
        # must be a development install.
        version = open(path).read().strip()
    else:
        # This is an end-user installation.
        if version_method == 'pkg_resources':
            version = pkg_resources.get_distribution('vespa-suite').version
        else:
            version = importlib.metadata.version('vespa-suite')

    return version


def get_windows_special_folder_path(folder_id):
    r"""Given a Windows CSIDL_XXX value (CSIDL = Constant Special ID List),
    returns the associated directory. For instance, passing 0x026
    (the value associated with CSIDL_PROGRAM_FILES) returns this:
       C:\Program Files\

    Some CSIDL constants are defined in constants.WindowsSpecialFolderIds.

    This function raises an error if called under a non-Windows OS.
    """
    # Getting Windows folders requires requires calling a function
    # provided by Windows itself. Thanks to Jay on stackoverflow for the
    # solution.
    # ref: http://stackoverflow.com/questions/626796/how-do-i-find-the-windows-common-application-data-folder-using-python

    # I import these here rather than at the top of the module because
    # I'm not sure windll & wintypes are available on non-Win platforms.
    import ctypes
    from ctypes import windll, wintypes

    _SHGetFolderPath = windll.shell32.SHGetFolderPathW
    _SHGetFolderPath.argtypes = [wintypes.HWND, ctypes.c_int,
                                 wintypes.HANDLE, wintypes.DWORD,
                                 wintypes.LPCWSTR]


    path_buffer = ctypes.create_unicode_buffer(wintypes.MAX_PATH)
    _SHGetFolderPath(0, folder_id, 0, 0, path_buffer)

    return path_buffer.value


def is_floatable(s):
    """True if the passed value can be turned into a float, False otherwise"""
    rc = False

    try:
        float(s)
        rc = True
    except ValueError:
        pass

    return rc


def is_intable(s):
    """True if the passed value can be turned into a int, False otherwise"""
    rc = False

    try:
        int(s)
        rc = True
    except ValueError:
        pass

    return rc


def is_iterable(an_object, include_strings=True):
    """ Returns True if an object is iterable, False otherwise. Iterable
    types include lists, tuples, dicts, strings, numpy arrays, etc.

    If include_strings is False, then strings are not considered iterable.
    This is useful because often we use is_iterable() to decide whether
    something is an object or a list of objects, and while one *can* loop
    over the contents of a string with 'for item in some_string', one
    doesn't typically want to.
    """
    try:
        iter(an_object)
        rc = True
    except TypeError:
        rc = False

    if rc and (not include_strings) and isinstance(an_object, str):
        rc = False

    return rc


def is_gzipped(f):
    """Given a file object or filename, returns True if the file appears to
    be in gzip format; False otherwise.

    The parameter can be a file object or a string indicating a filename.

    If the parameter is a file object, it must be opened with the mode 'rb'.
    It will still be open when the function returns, but the file pointer
    may be changed.
    """
    # Gzip files start with two bytes: ID1 and ID2.
    # "These have the fixed values ID1 = 31 (0x1f, \037),
    # ID2 = 139 (0x8b, \213), to identify the file as being in gzip format."
    #
    # ref: http://www.gzip.org/zlib/rfc-gzip.html#file-format
    GZIP_BYTES = '\x1f\x8b'

    please_close = False
    if hasattr(f, "read"):
        # It's already file-ish
        f.seek(0)
    else:
        # It's string-y
        f = open(f, "rb")
        please_close = True

    s = f.read(2)
    if please_close:
        f.close()

    return (s == GZIP_BYTES)


def iter_flatten(*args):
    """ Returns a flattened copy of an object if an object is iterable,
    Iterable types include lists, tuples, dicts, strings, numpy arrays, etc.

    This is useful for iterating through an entire list of lists one item
    at a time, but without losing the structure of the original list of lists.
    """
    for x in args:
        if hasattr(x, '__iter__'):
            for y in iter_flatten(*x):
                yield y
        else:
            yield x


def normalize_isotope_name(name):
    """Given a string that is an isotope name in some form (see examples
    below), attempts to normalize to the way we likes 'em, which is with the
    number first and the letter(s) capitalized.

    Returns the normalized string if successful, None otherwise.

    The function is not case- or whitespace-sensitive.

    All of these are valid input: "1H", "1h", "h1", " h2 ", " Li 7  ",
    "99Q", "00ZZ". The last two demonstrate that this function merely cares
    if the name is in the correct form; it does not check that the name
    exists in the database or that the atom exists in the periodic table.
    """
    number = ""
    letter = ""
    # Strip all whitespace (leading, trailing and inner).
    name = name.replace(" ", "").replace("\t", "").replace("\n", "")
    match = _ISOTOPE_REGEX.match(name)
    if match:
        # match.groups() will always return four groups, two of which will
        # be empty. After discarding the empties, the remaining two matches
        # are either (number, letter) or (letter, number).
        number, letter = [group for group in match.groups() if group]
        if number[0] not in _ALL_DIGITS:
            # oops, I have them reversed.
            number, letter = letter, number

        return number + letter.upper()
    else:
        return None


def normalize_newlines(a_string, target_newline="\n"):
    """ Given a string (can be Unicode or not), turns all newlines into
    the target newline which defaults to \n.
    """
    return _EOL_REGEX.sub(target_newline, a_string)


def safe_attribute_compare(this, that, attribute_name):
    """Compares the named attribute on the objects this and that. The
    attribute need not exist on either object.

    Returns True if --
       (a) Neither object has the attribute, OR
       (b) Both objects have the attribute and the value of the attr is equal.
    Returns False if --
       (c) One object has the attribute and the other doesn't, OR
       (d) Both objects have the attribute and the value of the attr differs.
    """
    this_has = hasattr(this, attribute_name)
    that_has = hasattr(that, attribute_name)

    # For starters, both objects must have the attr.
    equal = (this_has == that_has)

    if equal and this_has:
        # Both this and that have the attr. This covers cases (b) and (d).
        equal = (getattr(this, attribute_name) == getattr(that, attribute_name))
    #else:
        # They're not equal or they're both False. In both cases, there's
        # nothing further to be done; equal already holds the correct value.
        # This is cases (a) and (c).

    return equal


def safe_attribute_set(source, target, attribute_name):
    """Given two objects and an attribute name, copies the value of the
    attribute from source to target.

    Both objects must have the attribute in question or this function
    raises AttributeError.

    It's insurance for when we're copying attributes from objects that have
    different types and therefore might have different attribute names.

    """
    source_has = hasattr(source, attribute_name)
    target_has = hasattr(target, attribute_name)

    if source_has and target_has:
        setattr(target, attribute_name, getattr(source, attribute_name))
    else:
        raise AttributeError("""Attribute "%s" must exist on both source and destination""" % attribute_name)

def uuid():
    """Returns a random UUID. This function should be the source of all
    UUIDs generated by Vespa.
    """
    return str(uuid_module.uuid4())



def _test():
    # Unit tests for is_iterable()
    assert(is_iterable([]))
    assert(is_iterable(""))
    assert(is_iterable({}))
    assert(is_iterable(()))
    assert(is_iterable("", False) == False)


    # Unit tests for safe_attribute_compare()
    class Anonymous(object):
        pass

    this = Anonymous()
    that = Anonymous()

    this.foo = 42

    assert(safe_attribute_compare(this, that, "bar"))
    assert(not safe_attribute_compare(this, that, "foo"))
    that.foo = 42
    assert(safe_attribute_compare(this, that, "foo"))
    that.foo = 43
    assert(not safe_attribute_compare(this, that, "foo"))




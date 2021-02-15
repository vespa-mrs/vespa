# Python modules

import os
import logging
import sys
import errno
import os.path
import shutil
import time

# 3rd party modules
import wx

# Our modules
import vespa.common.configobj as configobj
import vespa.common.util.config as util_config
import vespa.common.util.logging_ as util_logging
import vespa.common.util.db_upgrader as db_upgrader
import vespa.common.util.misc as util_misc
import vespa.common.util.dir_list as dir_list
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.default_ini_file_content as default_ini_file_content
import vespa.common.exception_handler as exception_handler

# _ONE_YEAR measures a year in seconds. It's a little sloppy because it doesn't
# account for leap years and that kind of thing, but it's accurate enough for
# how we use it.
_ONE_YEAR = 365 * 24 * 60 * 60


def init_app(app_name):
    """
    The function that all Vespa apps must call to create their wx.App
    object (which this function returns). The caller must pass the app name.
    The app name will be used in messages displayed to the user by this code,
    so it should be something user-friendly like RFPulse, Simulation or
    Analysis.

    There's nothing magic about this code, but it does some important
    housekeeping that we don't want to duplicate in each app.

    Specifically, it does the following --
    - Creates the wx.App object
    - Ensures only one copy of this application is running. It displays a
      message and quits if the app is already running.
    - Sets the wx app name
    - Creates the .vespa object that hangs off of the wx.App object
    - Ensures that Vespa data dir exists. If not it creates it and a new
      database as well.
    - Ensures that Vespa log dir exists, creates it if not.
    - (Re)creates any of our INI files that might be missing.
    - Optionally sets a custom exception hook.
    - If the database is missing, gives the user the option of recreating it
      or exiting the app.

    """
    first_run = False

    # When running under Windows w/pythonw.exe, we have to expicitly send
    # stdout and stderr to the bit bucket.
    # See http://scion.duhs.duke.edu/vespa/project/ticket/44
    # Code below adapted from http://bugs.python.org/issue706263#msg97442
    if (util_misc.get_platform() == "windows") and \
       ("pythonw" in sys.executable.lower()):
       blackhole = open(os.devnull, 'w')
       sys.stdout = sys.stderr = blackhole

    app = wx.App(False)

    app.SetAppName(app_name)

    if util_misc.get_platform() == "windows":
        # Under Windows, all versions of Vespa prior to 0.4.0
        # created the data dir in a suboptimal location. This code
        # checks for such a directory and, if it finds one, copies it to the
        # new, preferred location. It also cleans up the old one if it is
        # not in use.
        # ref: http://scion.duhs.duke.edu/vespa/project/ticket/39
        roaming_dir = util_misc.get_windows_special_folder_path(util_misc.WindowsSpecialFolderIds.CSIDL_APPDATA)
        # Normally hardcoding strings is a bad idea, but here it is
        # appropriate. All versions of Vespa that used the roaming dir
        # used the name "Vespa" for the directory name.
        roaming_dir = os.path.join(roaming_dir, "Vespa")
        if os.path.exists(roaming_dir):
            # OK, the non-preferred roaming directory is still around. Let's
            # see about the preferred location.
            data_dir = util_misc.get_data_dir()
            if os.path.exists(data_dir):
                # The preferred location exists, which means that the directory
                # in the roaming dir is no longer used. Remove it.
                shutil.rmtree(roaming_dir, True)
            else:
                # The preferred location does *not* exist, so we leave the
                # roaming version alone for the moment (as a backup) and make
                # a copy of it in the preferred location.
                shutil.copytree(roaming_dir, data_dir)



    # Create the Vespa data directory if necessary
    data_dir = util_misc.get_data_dir()
    if not os.path.exists(data_dir):
        # This looks like the first time a Vespa app has been run on this
        # computer.
        first_run = True
        os.mkdir(data_dir)

    # Ensure there isn't another instance of this app running. This code
    # must come after the creation of the app object because if the except
    # clause is triggered it will attempt to show a message box, and it's
    # not valid to do that before the app object exists.
    # Also, it must come after the creation of the data dir because it
    # writes a lock file there.
    try:
        single_instance_enforcer = __SingleAppInstanceEnforcer(app_name)
    except _AlreadyRunningException:
        common_dialogs.message("%s is already running." % app_name)
        sys.exit(-1)

    # Create a directory for log files if necessary
    log_dir = os.path.join(data_dir, "logs")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    _clean_log_files(log_dir)

    # Vespa <= 0.5.3 used INI file section and element names that were
    # somewhat different from what we use now.
    # Translating one to the other invites madness since there are many small
    # differences. Instead, we just delete the old INI files and create new
    # ones.
    # To ease the pain of losing one's preferences, we save the paths out of
    # the existing INI files (if any) and reinsert them in the new INI files.
    saved_ini_paths = { }

    for ini_name in ("simulation", "analysis", "datasim", "pulse"):
        filename = os.path.join(data_dir, "%s.ini" % ini_name)
        if os.path.exists(filename):
            # The INI file exists; see if it is from vespa <= 0.5.3.
            content = open(filename, "rb").read()

            # I detect old-style INI files by looking for a comment string
            # that has been in the default INI content for as long as we've
            # had default_ini_file_content.py (2 years). It's unlikely to get
            # used again because we don't have a single foreground color
            # anymore.
            if "#    foreground_color = black" in str(content):
                # This is an old-school INI file.
                content = util_misc.normalize_newlines(content)
                content = content.split('\n')
                config = configobj.ConfigObj(content, encoding="utf-8")

                if "paths" in config:
                    saved_ini_paths[ini_name] = config["paths"]

                os.remove(filename)
            # else:
                # Nothing to do -- this is a modern INI file
        # else:
            # Nothing to do -- the file doesn't exist

    # (Re)create any missing INI files.
    for ini_name in ("vespa", "simulation", "analysis", "pulse", "datasim",
                     "analysis_import_menu_additions", ):
        filename = os.path.join(data_dir, "%s.ini" % ini_name)
        if not os.path.exists(filename):
            _create_ini_file(filename, ini_name)

            if ini_name in saved_ini_paths:
                # Restore these paths
                config = configobj.ConfigObj(filename, encoding="utf-8")
                config["paths"] = saved_ini_paths[ini_name]
                config.write()


    # We add to the app object an object of our own. The object doesn't do
    # anything other than create a place for us to attach app-level (global)
    # things. We could attach them directly to the app object but then we
    # run the risk of trashing some wx-specific attribute that has the same
    # name.
    class Anonymous(object):
        """A generic object to which one can attach arbitrary attrs"""
        pass

    app.vespa = Anonymous()

    # The singleton enforcer needs to survive for the life of the app so
    # we stash a reference to it in the app object.
    app.vespa.single_instance_enforcer = single_instance_enforcer

    # The database name is in vespa.ini
    config = util_config.VespaConfig()

    # Here we deal with the [debug] section of vespa.ini. See here --
    # http://scion.duhs.duke.edu/vespa/project/wiki/IniFiles

    # We set our custom exception hook as early as possible. Note that it
    # depends on the directory for log files and the existence of vespa.ini,
    # so we can't set the exception hook before those exist.
    # The exception hook is on by default but can be turned off by a flag in
    # vespa.ini.
    hook_exceptions = True
    if ("debug" in config) and ("hook_exceptions" in config["debug"]):
        hook_exceptions = config["debug"].as_bool("hook_exceptions")

    if hook_exceptions:
        sys.excepthook = exception_handler.exception_hook

    # numpy's response to IEEE 754 floating point errors is configurable.
    # On many platforms these errors are sent to the bit bucket. As
    # developers, we want a chance to hear them loud & clear.
    # See here for details:
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.seterr.html
    if ("debug" in config) and ("numpy_error_response" in config["debug"]):
        numpy_error_response = config["debug"]["numpy_error_response"]
        if numpy_error_response:
            import numpy
            numpy.seterr(all=numpy_error_response)

    db_path = os.path.join(data_dir, config["database"]["filename"])

    create_database = False

    if first_run:
        # We create the database unconditionally
        create_database = True
    else:
        # Ensure the database exists
        if not os.path.exists(db_path):
            # Uh-oh
            msg = "Sorry, %s can't find its database. It should be here:\n" % app_name
            msg += db_path
            msg += "\n\n%s can create a new database and continue, or " % app_name
            msg += "it can quit now. Would you like to create a new database?"
            if wx.YES == common_dialogs.message(msg,
                                                "%s Initialization Error" % app_name,
                                                common_dialogs.Q_YES_NO):
                create_database = True
            else:
                # No database? We can't continue.
                sys.exit(-1)


    if create_database:
        wx.BeginBusyCursor()
        logger = logging.getLogger(util_logging.Log.DATABASE)
        # Set up a handler for this log
        logger.addHandler(util_logging.DatabaseFileHandler())

        # I import create_database here instead of unconditionally
        # because it imports lots of modules that we won't otherwise
        # need.
        import vespa.common.create_database as create_database

        create_database.create_database(db_path, logger)

        wx.EndBusyCursor()

    # Upgrade the database if necessary
    db = db_upgrader.DatabaseUpgrader(db_path)
    db.upgrade()

    return app, db_path


class _AlreadyRunningException(Exception):
    """Indicates that this app is already running. Raised only by the
    __SingleAppInstanceEnforcer class.
    """
    def __init__(self):
        Exception(self)


class __SingleAppInstanceEnforcer(object):
    """A class that ensures only a single instance of the app named in
    app_name (a string passed to the contructor) can be created.

    I swiped most of the code from the URL below and modified it a little
    for Vespa.
    http://stackoverflow.com/questions/380870/python-single-instance-of-program
    """
    def __init__(self, app_name):
        already_running = False

        self.lockfile = "%s.lock" % app_name

        self.lockfile = os.path.join(util_misc.get_data_dir(), self.lockfile)

        if sys.platform.startswith("win"):
            # Windows
            try:
                # The file already exists. This may be because the previous
                # execution was interrupted, so we try to remove the file.
                if os.path.exists(self.lockfile):
                    os.unlink(self.lockfile)
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
            except OSError as e:
                if e.errno == errno.EACCES:
                    already_running = True
                else:
                    print(e.errno)
                    raise
        else:
            # non Windows. Note that fcntl isn't available on Windows
            import fcntl

            self.fp = open(self.lockfile, 'w')
            try:
                fcntl.lockf(self.fp, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except IOError:
                already_running = True

        if already_running:
            raise _AlreadyRunningException


    def __del__(self):
        import sys
        if sys.platform.startswith("win"):
            if hasattr(self, 'fd'):
                os.close(self.fd)
                os.unlink(self.lockfile)



def _clean_log_files(log_file_path):
    # This removes any log files that are > 1 year old. This ensures that our
    # log files don't grow infinitely.
    # It assumes the log file path exists, so don't call it until after that
    # path has been created.
    lister = dir_list.DirList(log_file_path, dir_list.FilterFile(),
                                             dir_list.FilterEndsWith(".txt"))

    # This temporary function returns True if a file is more than a year old,
    # False otherwise.
    now = time.time()
    f = lambda filename: (now - int(os.path.getctime(filename))) > _ONE_YEAR

    # Find all files that match
    filenames = [filename for filename in lister.fq_files if f(filename)]

    # Whack 'em.
    list(map(os.remove, filenames))


def _create_ini_file(filename, ini_name):
    # Creates the INI file for an app (rfpulse, sim, or analysis) using
    # default content.
    content = default_ini_file_content.DEFAULT_INI_FILE_CONTENT[ini_name]

    open(filename, "w").write(content)



# Python modules

import sys
import os
import stat
import subprocess
import tempfile

# Our modules
import vespa.common.util.misc as util_misc

"""Creates double-clickable shortcuts for all Vespa apps.
It will only work when Vespa has been successfully installed.

This code is meant to be run like so:
python -m 'vespa.create_shortcuts'
"""


# APPLICATIONS lists the directories that contain the Vespa apps installed
# by setup.py.
APPLICATIONS = ("pulse", "simulation", "analysis", "datasim")

# Vespa supports OS X, Windows and Linux, and all three have fairly different
# ways of creating executable shortcuts.

platform = sys.platform.lower()
if "linux" in platform:
    platform = "linux"
elif "darwin" in platform:
    platform = "osx"
elif "win32" in platform:
    platform = "windows"

vespa_install_path = util_misc.get_vespa_install_directory()
python_path = sys.executable

done_msg = "Done! Vespa shortcuts have been created on your desktop."

if platform == "osx":
    # Under OS X we create shortcuts in the standard Applications folder. The shortcuts are
    # directories with the .app extension which makes them look and behave like resgular OS X
    # apps. They are bare bones apps, though, because they only contain a single 2-line shell
    # script.
    target_path = '/Applications/Vespa'
    done_msg = "Done! Vespa shortcuts have been created in your Applications folder."

    if not os.path.exists(target_path):
        os.mkdir(target_path)

    for application in APPLICATIONS:
        capitalized_name = application.capitalize()
        print("Creating a shortcut for %s..." % capitalized_name)

        # Create directory tree, e.g. /Applications/Vespa/Simulation.app/Contents/MacOS
        app_path = os.path.join(target_path, capitalized_name + ".app", "Contents", "MacOS")
        if not os.path.exists(app_path):
            os.makedirs(app_path)
        app_path = os.path.join(app_path, capitalized_name)

        # Create a command string that will run this app. The path to Python must be
        # fully qualified.
        path = os.path.join(vespa_install_path, application, "src", "main.py")
        executable = sys.executable
        if 'conda' in executable:
            # Mini/Anaconda on the Mac requires pythonw to launch GUI apps.
            executable += 'w'
        command = '"%s" "%s"' % (executable, path)

        with open(app_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(command)

        os.chmod(app_path, 0o755)

elif platform == "linux":
    # Get some environment vars.
    home = os.getenv("HOME")

    for application in APPLICATIONS:
        print("Creating a shortcut for %s..." % application.capitalize())

        # Create a command string that will run this app
        app_path = os.path.join(vespa_install_path, application, "src", "main.py")

        app_name = application.capitalize()

        # Gnome & KDE both support .desktop files. They're plain text
        # files with a well-documented standard. What a crazy idea!
        # http://standards.freedesktop.org/desktop-entry-spec/desktop-entry-spec-latest.html

        # Note that under Ubuntu (and perhaps other distros), pythonw is
        # not available by default. Under OS X and Windows, it is.
        content = """
#!/usr/bin/env xdg-open

[Desktop Entry]
Encoding=UTF-8
Version=1.0
Type=Application
Terminal=false
Icon[en_US]={python_path}
Name[en_US]={app_name}
Exec={python_path} {app_path}
Name={app_name}
Icon=python
""".format(python_path=python_path, app_name=app_name, app_path=app_path)

        path = home + ("/Desktop/%s.desktop" % app_name)
        open(path, "w").write(content)

        # Make the file executable.
        mode = os.stat(path)[0]
        mode |= stat.S_IXUSR
        os.chmod(path, mode)

elif platform == "windows":
    # Under Windows, executable shortcuts are .lnk files which are an
    # undocumented binary format. Thanks, Microsoft! The only safe way
    # to create them is to call functions in Windows DLLs. The easiest (?)
    # way for us to call those functions is via VBScript.
    # We write our VBScript to a temp file and then ask Windows to execute
    # it. Hopefully executing a .vbs won't cause any Windows security
    # thingies to freak out.

    # Get the path to the Python executable and change it to the pythonw
    # EXE if possible. In contrast to python.exe, pythonw.exe doesn't
    # open a console window before launching the app, so it looks a bit
    # nicer. The downside is that anything printed to stdout or stderr
    # will go into the bit bucket.
    # Under non-Windows platforms, there's no difference between python
    # and pythonw.
    if python_path.endswith("python.exe"):
        python_path = python_path[:-len("python.exe")] + "pythonw.exe"
    # else:
        # This is very unexpected, better not mess with it.

    for application in APPLICATIONS:
        print("Creating a shortcut for %s..." % application.capitalize())

        # Create the path to main.py for this app
        app_path = os.path.join(vespa_install_path, application, "src", "main.py")

        name = application.capitalize()

        # This is the VBScript.
        content = """
WScript.Quit Main

Function Main
Set shell = CreateObject("WScript.Shell")

With shell.CreateShortcut(shell.SpecialFolders("Desktop") & "\\%s.lnk")
.TargetPath = chr(34) + "%s" + chr(34)
.Arguments = chr(34) + "%s" + chr(34)
.WindowStyle = 1
.Save
End With
End Function
""" % (name, python_path, app_path)

        # I create a temp file from which to execute the script. The temp
        # file must have a .vbs extension.
        fd, filename = tempfile.mkstemp(".vbs")
        os.write(fd, content)
        os.close(fd)

        # Execute it
        os.system(filename)

        # Clean up
        os.remove(filename)

print(done_msg)

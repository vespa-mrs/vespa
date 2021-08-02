

import os
import sys
import collections
import tempfile
import subprocess   

# This script is meant to be run before Vespa is installed, so you must *not* import any
# Vespa modules.

"""
This script has 3 similar modes, all of which are based on its ability to sniff out and examine
all the Pythons present in the directories present in the PATH environment variable.

The 3 modes are --
   1. As a Vespa pre-install step, this script can report to the user which Pythons are
   Vespa-compatible and allow the user to select one. If that Python is not already first in the
   PATH, it prints a shell command that will temporarily make it first in the PATH.

   2. As a Vespa pre-upgrade step, this script can report to the user the Python where Vespa is
   installed. If that Python is not already first in the PATH, this script prints a shell
   command that will temporarily make it first in the PATH.

   3. Print a list of installed Pythons, including which ones for which Vespa is installed.

One can pass a mode number (1, 2, or 3) as an argument when invoking the script.
If no argument is passed, the script prompts for a selection.

This script recognizes 32-bit Python, but not 32-bit Vespa. In other words, if 32-bit Vespa is
installed for a Python, that Vespa install will be invisible to this script.

This script knows about the OS X system Python and discourages users from using it.

This script is and must remain compatible with Python 2 and 3. It will typically be executed by
the first Python in a user's PATH, and there's no guarantee that will be Python 2.

It must also be runnable from anywhere on the file system (e.g. a user's Downloads folder).
"""

try:
    raw_input
except NameError:
    # This is Python 3. Make raw_input an alias for Python 3's input() builtin.
    raw_input = input


# PATH_SEPARATOR is the delimiter in the PATH environment variable
PATH_SEPARATOR = ';' if (sys.platform == 'win32') else ':'


class Python(object):
    """Represents a Python executable.

    - path is the absolute, normalized path to a Python executable, including the executable
    itself (e.g. C:\Python27\python.exe)
    - bit_mode is an int, either 32 or 64
    - version is the 3-digit version string, e.g. '2.7.11'
    - is_vespa_installed is True if Vespa is installed for this Python, False otherwise.
    It cannot be True when is_vespa_compatible is False.
    """
    def __init__(self, path, bit_mode, version):
        self.path = path
        self.bit_mode = bit_mode
        self.version = version
        # is_vespa_installed is set after init.
        self.is_vespa_installed = False

    @property
    def is_osx_system_python(self):
        return (sys.platform == 'darwin') and \
               (self.path.startswith('/System/Library/Frameworks/Python.framework') or
                self.path.startswith('/usr/bin/'))

    @property
    def is_vespa_compatible(self):
        return (self.bit_mode == 64) and self.version.startswith('2.7')

    def __str__(self):
        return "{}-bit Python {} at {}".format(self.bit_mode, self.version, self.path)


def is_powershell():
    """True if the default shell is Windows Powershell, False otherwise."""
    if sys.platform == 'win32':
        # Both DOS and Powershell support the echo command. Under DOS, %foo% is special syntax
        # that displays the value of the environment variable foo. Under Powershell it has no
        # special meaning and so it just returns the same string. Therefore, if 'echo %PATH%'
        # returns '%PATH%', the shell is Powershell, otherwise the shell is DOS.
        return '%PATH%' == subprocess.check_output(('echo', '%PATH%'), shell=True)
    else:
        # Not Windows
        return False


def build_path_command(path_to_prepend):
    """Given a path, returns a string that, executed at the OS command prompt, will prepend the
    path to the PATH environment variable.
    """
    path_to_prepend += PATH_SEPARATOR
    if sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        # This works in bash, not sure about other shells. bash is the default under OS X and
        # Ubuntu, and the latter is the only Linux we officially support.
        cmd = 'export "PATH={}$PATH"'.format(path_to_prepend)
    else:
        # Windows. This could be cmd.exe (DOS) or Powershell.
        if is_powershell():
            cmd = '$env:PATH = "{}$env:PATH'.format(path_to_prepend)
        else:
            cmd = 'PATH = "{}%PATH%"'.format(path_to_prepend)

    return cmd


def prettify_numbers(numbers, conjunction='and'):
    """Given a list of 0 or more ints, returns a string that presents them nicely, e.g.
    [1, 3, 5] ==> '1, 3, and 5'
    [1, 5] ==> '1 and 5'
    [2] ==> '2'
    [] ==> ''

    The conjunction is present only if there's multiple numbers ('1, 3, and 5' vs. '1, 3, or 5').
    """
    pretty = ''

    if len(numbers) == 0:
        # Nothing to do
        pass
    elif len(numbers) == 1:
        pretty += str(numbers[0])
    elif len(numbers) == 2:
        pretty += '{} {} {}'.format(numbers[0], conjunction, numbers[1])
    else:
        # len(numbers) > 2
        pretty += ', '.join([str(x) for x in numbers[:-1]])
        pretty += ", {} {}".format(conjunction, str(numbers[-1]))

    return pretty


def normalize_path(path):
    """Given an absolute path, normalizes it in two ways and returns the normalized path.

    The first normalization is that symlinks are followed until they return a real path.

    The second normalization is that the path is run through os.path.normpath().
    """
    while os.path.islink(path):
        # This is a symlink (only happens on *nix). Symlinks can be relative, and if they are,
        # they're relative to the directory in which they reside. e.g. /usr/local/bin/python
        # might be a symlink defined as ../../bin/python which resolves to /usr/bin/python.
        # The two lines below make sense of relative (and absolute) symlinks.
        dirname = os.path.dirname(path)
        path = os.path.join(dirname, os.readlink(path))

    # Path is absolute but possibly messy; e.g. /usr/local/bin/../../bin/python.
    # normpath() cleans that up.
    return os.path.normpath(path)


def find_pythons():
    """Returns a list describing the Pythons installed on this system.

    Each Python installation is represented by an instance of the Python class defined in
    this module.

    This only looks for Python in directories in the PATH.
    """
    # Some Python 3 installations name the executable 'python3[.exe]', not just 'python[.exe]'.
    if sys.platform == 'win32':
        executable_names = ('python.exe', 'python3.exe')
    else:
        executable_names = ('python', 'python3')

    # I store the Pythons I find as a dict so I can easily filter duplicates.
    pythons = collections.OrderedDict()
    path = os.environ['PATH'].split(PATH_SEPARATOR)
    for directory in path:
        for executable_name in executable_names:
            possible_python = os.path.join(directory, executable_name)
            if os.path.exists(possible_python):
                # Under *nix, one Python can exist in multiple locations in the file system due
                # to symlinks. e.g. on my Mac, the Python 3.5 executable is here:
                #    /Library/Frameworks/Python.framework/Versions/3.5/bin/python3
                # It's also pointed to by this symlink:
                #    /usr/local/bin/python3
                # Both of these directories appear in my PATH. So as to minimize confusion when
                # reporting Pythons to the user, it's important that this code recognizes that the
                # two entries actually represent just one Python. normalize_path() helps with that.
                possible_python = normalize_path(possible_python)

                if possible_python in pythons:
                    # Nevermind, I've seen this one already.
                    pass
                else:
                    # Try to execute it
                    args = (possible_python, '-c', "import struct; print(8 * struct.calcsize('P'))")
                    try:
                        bit_mode = subprocess.check_output(args)
                    except OSError:
                        # This is a bit unexpected, but not completely so we ignore it.

                        # Note that you will also get this exception if your command contains a
                        # Python syntax error (e.g. unbalanced parens, or Python 2-specific syntax
                        # executed under Python 3), so if you're doing development it's usually a
                        # good idea to add a 'raise' statement here to re-raise errors so you don't
                        # obscure an ordinary coding error.
                        bit_mode = None

                    if bit_mode:
                        bit_mode = int(bit_mode)

                        args = (possible_python, '-c',
                                "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))")
                        try:
                            version = subprocess.check_output(args)
                        except OSError:
                            # See comment above about handling OSError.
                            version = None

                    if bit_mode and version:
                        # We found a Python. Save it.
                        # subprocess.check_output() returns a byte string. Under Python 3 it's
                        # important to decode it to a real (Unicode) string.
                        version = version.decode('ASCII').strip()
                        python = Python(possible_python, bit_mode, version)
                        pythons[possible_python] = python

                        if python.is_vespa_compatible:
                            # See if Vespa is installed under this Python.
                            args = (python.path, '-c', 'import vespa')
                            try:
                                # Some non-obvious stuff here. I call check_output() and redirect
                                # stderr to stdout not because I care about the output (whether or
                                # not an exception is raised tells me all I need to know) but
                                # because if I don't capture the output, it will be displayed to the
                                # user and it looks bad.
                                # Also, I have to set the CWD because Python adds the CWD to
                                # sys.path when it starts. As a developer, I often invoke this
                                # script when my CWD is something where 'import vespa' succeeds by
                                # importing it from the CWD (rather than site-packages), giving me
                                # a false positive.
                                subprocess.check_output(args, cwd=tempfile.gettempdir(),
                                                        stderr=subprocess.STDOUT)
                            except subprocess.CalledProcessError:
                                python.is_vespa_installed = False
                            else:
                                python.is_vespa_installed = True

    return list(pythons.values())


def print_install_compatibility(pythons):
    """Given a non-empty list of Pythons, implements mode 1 (see file docstring)"""
    python = pythons[0]

    if python.is_vespa_compatible and not python.is_osx_system_python:
        # The first Python in the PATH is Vespa-compatible, so all is well.
        print("Everything looks good! You're ready to proceed with the Vespa installation.")
    else:
        compatible = []
        system_pythons = []
        msg = "I found the following Pythons installed on your system:\n"
        for i, python in enumerate(pythons, 1):
            if python.is_vespa_compatible:
                compatible.append(i)

                if python.is_osx_system_python:
                    system_pythons.append(i)
            # else:
                # If the Python isn't Vespa-compatible, we don't care if it's the system Python
                # because the user can't choose it anyway.

            msg += "   {}: {}\n".format(i, python)

        if compatible:
            msg += "\nVespa can use Python {}. ".format(prettify_numbers(compatible))

            if system_pythons:
                msg += "\nHowever, Python {} ".format(prettify_numbers(system_pythons))
                msg += "is/are the OS X system Python. "
                msg += "We recommend you do NOT use the system Python for Vespa.\n\n"

            if (1 in compatible) and (not pythons[0].is_osx_system_python):
                # The first Python in the PATH is Vespa-compatible and it's not the OS X system
                # Python. No need to bug the user with any questions.
                print("Everything looks good! You're ready to proceed with the Vespa installation.")
            else:
                msg += "Which Python do you want to use for Vespa?\n"
                print(msg)

                compatible = [str(x) for x in compatible]
                prompt = "Enter {}, or x to exit: ".format(prettify_numbers(compatible, 'or'))

                selection = input(prompt)
                while (selection != 'x') and (selection not in compatible):
                    print("Sorry, I didn't understand that.")
                    selection = input(prompt)

                if selection == 'x':
                    # OK, they want to quit. Say something nice.
                    msg = "Please ensure you're using an appropriate Python before you proceed "
                    msg += "with the Vespa installation."
                    print(msg)
                    sys.exit(0)
                else:
                    selection = int(selection)

                    python = pythons[selection - 1]

                    msg = "\nPlease copy and paste this command to temporarily change your PATH.\n"
                    msg += "Next, proceed with the Vespa install. Be sure to use this command \n"
                    msg += "prompt for all of the installation steps."
                    print(msg)

                    print(('\n' + build_path_command(os.path.dirname(python.path)) + '\n'))
        else:
            # No Vespa-compatible Pythons found.
            print("I can't find a Vespa-compatible Python on this computer.")
            print('Please install a Vespa-compatible Python before you continue.')
            print('See the Vespa Web site for more information:')
            print('https://scion.duhs.duke.edu/vespa/project')


def print_upgrade_compatibility(pythons):
    """Given a non-empty list of Pythons, implements mode 2 (see file docstring).

    This assumes that Vespa is installed in only one location, and it ignores 32-bit Vespa.
    """
    # Find the Python with Vespa installed. (There might not be one at all.)
    python = pythons[0]

    if python.is_vespa_installed:
        # Vespa is installed under the first Python in the PATH, so all is well.
        print("Everything looks good! You're ready to proceed with the Vespa upgrade.")
    else:
        pythons = [python for python in pythons if python.is_vespa_installed]

        if pythons:
            python = pythons[0]

            msg = "Please copy and paste this command to temporarily change your PATH.\n"
            msg += "Next, proceed with the Vespa upgrade. Be sure to use this command prompt\n"
            msg += "for the upgrade step."
            print(msg)

            print(('\n' + build_path_command(os.path.dirname(python.path)) + '\n'))
        else:
            # Yikes; user is asking to upgrade Vespa but I can't find it at all.
            print("I can't find Vespa installed on this computer.")


def print_all_pythons(pythons):
    """Given a non-empty list of Pythons, implements mode 3 (see file docstring)."""
    print("I found the following Pythons installed on your system:")
    vespa_installations = []
    for i, python in enumerate(pythons, 1):
        print(("   {}: {}".format(i, python)))
        if python.is_vespa_installed:
            vespa_installations.append(i)

    print('')
    if vespa_installations:
        print(("Vespa is installed for Python {}.".format(prettify_numbers(vespa_installations))))
    else:
        print("Vespa is not installed for any of these Pythons.")


# +++++++++++++++++++++     main starts here     +++++++++++++++++++++

options = ['1', '2', '3']

if (len(sys.argv) > 1) and (sys.argv[1] in options):
    selection = int(sys.argv[1])
else:
    msg = "What would you like to do?\n"
    msg += "1) Check that your Python is Vespa-compatible\n"
    msg += "2) Check that your Python is ready for a Vespa upgrade\n"
    msg += "3) See a list of the Pythons installed on this computer\n"
    print(msg)
    prompt = "Enter {}, or x to exit: ".format(prettify_numbers(options, 'or'))
    selection = input(prompt)
    while (selection != 'x') and (selection not in options):
        print("Sorry, I didn't understand that.")
        selection = input(prompt)

    if selection == 'x':
        sys.exit(0)
    else:
        selection = int(selection)

pythons = find_pythons()

if not pythons:
    # We found 0 Pythons! This is very unexpected. If there's no Python in the PATH, then how
    # was this script invoked? Either the user invoked it with a full path to Python, or they
    # executed it from within the Python directory. Both of these cases are out of the scope of
    # what this script is expected to handle.
    print("I didn't find any Pythons installed on your system.")
    sys.exit(0)
else:
    if selection == 1:
        print_install_compatibility(pythons)
    elif selection == 2:
        print_upgrade_compatibility(pythons)
    else:
        print_all_pythons(pythons)

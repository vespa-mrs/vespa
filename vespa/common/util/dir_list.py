# Python modules

import fnmatch
import os
import itertools


class DirList(object):
    """A convenience class that gets a directory listing and provides
    various views of that listing. The listing can be filtered in a number
    of ways.

    At its simplest, a DirList instance is unfiltered and contains everything
    in a path as returned by os.listdir(). The following attributes/properties
    are available; all return a sorted (possibly empty) list --

    names - The names of all files and directories
    files - The names of the files; directories are not included
    dirs  - The names of the directories; files are not included
    fq_names, fq_files, fq_dirs -
            Same as above, but fully-qualified, meaning that the path
            passed to the constructor is prefixed to each name.

    In addition, the method fq(name) returns a fully-qualified version of
    name by prefixing the path given in the constructor.

    The DirList constructor accepts zero or more instances of the FilterXxx
    classes defined and described below. The filters are ANDed together. 
    
    Some examples --
    
    - To list all the directories in your music folder:
    DirList("./Music", FilterDir())
    
    - To list all the JPEGs in your artwork folder:
    DirList("./my_fabulous_art", FilterFnMatch("*.jpg"))
    
    - To fetch everything but the .pyc files in the "source" directory:
    DirList("./source", FilterEndsWith(".pyc", False))

    - To fetch files starting with "a" or "A" in the current directory:
    DirList(".", FilterFile(), FilterFnMatch("a*"))
    
    - To fetch directories starting with a number in the current directory:
    DirList(".", FilterDir(), FilterRegex(re.compile("^[\d]")))
    """

    def __init__(self, path, *args):
        self.path = path

        if not os.path.exists(path):
            raise ValueError("Path does not exist: %s" % path)

        self.names = sorted(os.listdir(path))

        # For each filter provided, I construct a function that compares a
        # filename to the filter's criterion and returns True or False.
        # I then apply this function to every name in self.names and thus
        # winnow the list.
        for filter_ in args:
            if isinstance(filter_, FilterStartsWith):
                filter_function = \
                    lambda name: name.startswith(filter_.criterion)
            elif isinstance(filter_, FilterEndsWith):
                filter_function = \
                    lambda name: name.endswith(filter_.criterion)
            elif isinstance(filter_, FilterFnMatch):
                filter_function = \
                    lambda name: fnmatch.fnmatch(name, filter_.criterion)
            elif isinstance(filter_, FilterFnMatchCase):
                filter_function = \
                    lambda name: fnmatch.fnmatchcase(name, filter_.criterion)
            elif isinstance(filter_, FilterRegex):
                filter_function = \
                    lambda name: filter_.criterion.search(name)
            elif isinstance(filter_, FilterFile):
                filter_function = \
                    lambda name: os.path.isfile(self.fq(name))
            elif isinstance(filter_, FilterDir):
                filter_function = \
                    lambda name: os.path.isdir(self.fq(name))
            else:
                filter_function = None

            if filter_function:
                if filter_.include_only:
                    self.names = list(filter(filter_function, self.names))
                else:
                    self.names = list(itertools.filterfalse(filter_function,
                                                             self.names))


    @property
    def files(self):
        """Returns files in the list."""
        f = lambda name: os.path.isfile(self.fq(name))
        return list(filter(f, self.names))


    @property
    def dirs(self):
        """Returns directories in the list."""
        f = lambda name: os.path.isdir(self.fq(name))
        return list(filter(f, self.names))


    @property
    def fq_names(self):
        """Returns a list of fully-qualified names."""
        return [self.fq(name) for name in self.names]


    @property
    def fq_files(self):
        """Returns a list of fully-qualified files from the list."""
        return list(filter(os.path.isfile, self.fq_names))


    @property
    def fq_dirs(self):
        """Returns a list of fully-qualified directories from the list."""
        return list(filter(os.path.isdir, self.fq_names))


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        return str(self.names)


    def fq(self, name):
        """Fully qualifies the name given by prepending the path that was
        passed in the constructor."""
        return os.path.join(self.path, name)



class Filter(object):
    """Base class for all the other filters. Don't use this as a filter on 
    it's own; it doesn't do anything. If you write your own filter, you should
    use this as the base class.

    Each filter accepts a criterion and a boolean include_only flag. The
    criterion differs from one filter to another and is explained in the
    filter's docstring.
    
    The include_only param is the same for each filter. When True (the
    default), only files matching the criterion are included in the DirList
    and all others are discarded. (This is similar to a whitelist.) When
    include_only=False, all files are included in the DirList except for those
    matching the criterion. This is similar to a blacklist.
    """
    def __init__(self, criterion, include_only=True):
        self.criterion = criterion
        self.include_only = include_only

class FilterStartsWith(Filter):
    """Filters if the name starts with the criterion.
    
    For example, to include only files starting with "rabbit" --
        FilterStartsWith("rabbit")
    To exclude all files starting with a dot --
        FilterStartsWith('.', False)
    """
    pass

class FilterEndsWith(Filter):
    """Filters if the name ends with the criterion.

    For example, to include text files --
        FilterEndsWith(".txt")
    To exclude compiled Python files --
        FilterEndsWith(".pyc", False)
    """
    pass

class FilterFnMatch(Filter):
    """Filters if the name matches according to the standard library's
    fnmatch.fnmatch(). 

    For example, to include PNG files --
        FilterFnMatch("*.png")
    To exclude all Python files --
        FilterFnMatch("*.py*", False)
    """
    pass

class FilterFnMatchCase(Filter):
    """Filters if the name matches according to the standard library's
    fnmatch.fnmatchcase(). Similar to FilterFnMatch.
    """
    pass

class FilterRegex(Filter):
    """Filters if the name matches the criterion which must be a compiled
    regular expression. Each name in the directory listing is compared
    using regex.search().
    """
    pass

class FilterFile(Filter):
    """Filters if the name refers to a file."""
    def __init__(self, include_only=True):
        Filter.__init__(self, None, include_only)

class FilterDir(Filter):
    """Filters if the name refers to a directory."""
    def __init__(self, include_only=True):
        Filter.__init__(self, None, include_only)




# FIXME - PS - would be nice to have some tests, but since this code relies on 
# what's reported by the filesystem, we'd have to create a temp directory,
# create a bunch of known files & directories in that temp directory, run
# tests against that and then clean up afterward. It's a lot to do right now,
# sorry!
